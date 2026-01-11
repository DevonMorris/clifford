# PRD-14.7: CLI Tool & Build Integration

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Provide a command-line tool for generating code and integrating with build systems

## Design Decision: CLI over Proc Macro

We chose a **CLI tool** over a proc macro for several reasons:

| Aspect | CLI | Proc Macro |
|--------|-----|------------|
| **Generation cost** | Once, then committed | Every build |
| **Build times** | Normal | Slower compiles |
| **Debuggability** | Can inspect generated code | Hidden in target/ |
| **IDE support** | Full (real .rs files) | Limited |
| **Incremental builds** | Works normally | Regenerates each time |
| **CI caching** | Effective | Less effective |
| **Error messages** | Clear file:line | Macro expansion errors |

The generated code is committed to the repository, meaning:
- Builds are fast (no generation step)
- IDE navigation and autocomplete work perfectly
- Code review shows exactly what will be compiled
- CI only regenerates when specs change

## CLI Interface

### Basic Usage

```bash
# Generate a single algebra
clifford-codegen generate algebras/euclidean3.toml -o src/generated/euclidean3/

# Generate multiple algebras
clifford-codegen generate algebras/*.toml -o src/generated/

# Verify generated code matches Multivector
clifford-codegen verify src/generated/euclidean3/

# Show what would be generated (dry run)
clifford-codegen generate algebras/pga3.toml --dry-run

# List all products for a type
clifford-codegen products algebras/cga3.toml --type Motor
```

### Full Command Reference

```
clifford-codegen

USAGE:
    clifford-codegen <SUBCOMMAND>

SUBCOMMANDS:
    generate    Generate Rust code from algebra specifications
    verify      Verify generated code matches Multivector
    products    List products for an algebra
    blades      List blades and their indices
    help        Print help information

---

clifford-codegen generate

USAGE:
    clifford-codegen generate [OPTIONS] <SPEC>...

ARGS:
    <SPEC>...    Path(s) to algebra specification files (.toml)

OPTIONS:
    -o, --output <DIR>        Output directory (default: src/generated/{algebra_name})
    --nalgebra                Include nalgebra conversions
    --serde                   Include serde derives
    --no-tests                Skip test generation
    --no-arbitrary            Skip Arbitrary implementations
    --dry-run                 Show what would be generated without writing files
    --force                   Overwrite existing files
    -v, --verbose             Verbose output
    -h, --help                Print help information

---

clifford-codegen verify

USAGE:
    clifford-codegen verify <PATH>

ARGS:
    <PATH>    Path to generated module directory

OPTIONS:
    --spec <FILE>    Original specification file (for re-verification)
    -h, --help       Print help information

---

clifford-codegen products

USAGE:
    clifford-codegen products [OPTIONS] <SPEC>

ARGS:
    <SPEC>    Path to algebra specification file

OPTIONS:
    --type <TYPE>     Show products for specific type
    --product <KIND>  Filter by product kind (geometric, outer, inner, sandwich)
    --table           Output as markdown table
    -h, --help        Print help information

---

clifford-codegen blades

USAGE:
    clifford-codegen blades <SPEC>

ARGS:
    <SPEC>    Path to algebra specification file

OPTIONS:
    --grade <N>     Show only blades of grade N
    --format <FMT>  Output format: table, json
    -h, --help      Print help information
```

## Implementation

### Main Entry Point

```rust
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "clifford-codegen")]
#[command(about = "Generate optimized geometric algebra code")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate Rust code from algebra specifications
    Generate {
        /// Specification file(s)
        #[arg(required = true)]
        specs: Vec<PathBuf>,

        /// Output directory
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Include nalgebra conversions
        #[arg(long)]
        nalgebra: bool,

        /// Include serde derives
        #[arg(long)]
        serde: bool,

        /// Skip test generation
        #[arg(long)]
        no_tests: bool,

        /// Skip Arbitrary implementations
        #[arg(long)]
        no_arbitrary: bool,

        /// Show what would be generated
        #[arg(long)]
        dry_run: bool,

        /// Overwrite existing files
        #[arg(long)]
        force: bool,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Verify generated code matches Multivector
    Verify {
        /// Path to generated module
        path: PathBuf,

        /// Original specification
        #[arg(long)]
        spec: Option<PathBuf>,
    },

    /// List products for an algebra
    Products {
        /// Specification file
        spec: PathBuf,

        /// Filter by type
        #[arg(long, name = "type")]
        type_name: Option<String>,

        /// Filter by product kind
        #[arg(long)]
        product: Option<String>,

        /// Output as markdown table
        #[arg(long)]
        table: bool,
    },

    /// List blades and indices
    Blades {
        /// Specification file
        spec: PathBuf,

        /// Filter by grade
        #[arg(long)]
        grade: Option<usize>,

        /// Output format
        #[arg(long, default_value = "table")]
        format: String,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Generate {
            specs,
            output,
            nalgebra,
            serde,
            no_tests,
            no_arbitrary,
            dry_run,
            force,
            verbose,
        } => {
            for spec_path in specs {
                generate(&spec_path, &GenerateOptions {
                    output: output.clone(),
                    nalgebra,
                    serde,
                    tests: !no_tests,
                    arbitrary: !no_arbitrary,
                    dry_run,
                    force,
                    verbose,
                })?;
            }
        }
        Commands::Verify { path, spec } => {
            verify(&path, spec.as_deref())?;
        }
        Commands::Products { spec, type_name, product, table } => {
            list_products(&spec, type_name.as_deref(), product.as_deref(), table)?;
        }
        Commands::Blades { spec, grade, format } => {
            list_blades(&spec, grade, &format)?;
        }
    }

    Ok(())
}
```

### Generate Command

```rust
struct GenerateOptions {
    output: Option<PathBuf>,
    nalgebra: bool,
    serde: bool,
    tests: bool,
    arbitrary: bool,
    dry_run: bool,
    force: bool,
    verbose: bool,
}

fn generate(spec_path: &Path, options: &GenerateOptions) -> Result<()> {
    // Parse specification
    let spec_content = std::fs::read_to_string(spec_path)?;
    let spec = parse_spec(&spec_content)?;

    if options.verbose {
        println!("Generating algebra: {}", spec.name);
        println!("  Signature: ({}, {}, {})", spec.signature.p, spec.signature.q, spec.signature.r);
        println!("  Types: {}", spec.types.len());
    }

    // Determine output directory
    let output_dir = options.output.clone()
        .unwrap_or_else(|| PathBuf::from(format!("src/generated/{}", spec.name)));

    // Check if output exists
    if output_dir.exists() && !options.force && !options.dry_run {
        return Err(anyhow!("Output directory exists. Use --force to overwrite."));
    }

    // Build algebra
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
    let table = ProductTable::new(&algebra);

    // Generate all files
    let files = vec![
        ("mod.rs", generate_mod(&spec)),
        ("types.rs", generate_types(&spec, &algebra)),
        ("constrained.rs", generate_constrained(&spec, &algebra)),
        ("products.rs", generate_products(&spec, &algebra, &table)),
        ("ops.rs", generate_ops(&spec)),
        ("approx.rs", generate_approx(&spec)),
        ("conversions.rs", generate_conversions(&spec, &algebra)),
    ];

    let optional_files = vec![
        (options.arbitrary, "arbitrary.rs", generate_arbitrary(&spec)),
        (options.tests, "tests.rs", generate_tests(&spec, &algebra)),
        (options.nalgebra, "nalgebra.rs", generate_nalgebra(&spec)),
    ];

    if options.dry_run {
        println!("Would generate:");
        for (name, _) in &files {
            println!("  {}/{}", output_dir.display(), name);
        }
        for (enabled, name, _) in &optional_files {
            if *enabled {
                println!("  {}/{}", output_dir.display(), name);
            }
        }
        return Ok(());
    }

    // Create output directory
    std::fs::create_dir_all(&output_dir)?;

    // Write files
    for (name, tokens) in files {
        let path = output_dir.join(name);
        write_formatted(&path, tokens)?;
        if options.verbose {
            println!("  Wrote: {}", path.display());
        }
    }

    for (enabled, name, tokens) in optional_files {
        if enabled {
            let path = output_dir.join(name);
            write_formatted(&path, tokens)?;
            if options.verbose {
                println!("  Wrote: {}", path.display());
            }
        }
    }

    println!("Generated {} in {}", spec.name, output_dir.display());

    Ok(())
}

fn write_formatted(path: &Path, tokens: TokenStream) -> Result<()> {
    // Convert to string
    let code = tokens.to_string();

    // Write to file
    std::fs::write(path, &code)?;

    // Run rustfmt
    std::process::Command::new("rustfmt")
        .arg(path)
        .status()?;

    Ok(())
}
```

### Verify Command

```rust
fn verify(path: &Path, spec_path: Option<&Path>) -> Result<()> {
    // Load the generated module
    let mod_rs = path.join("mod.rs");
    if !mod_rs.exists() {
        return Err(anyhow!("Not a generated module: {}", path.display()));
    }

    // If spec provided, regenerate and compare
    if let Some(spec_path) = spec_path {
        let spec = parse_spec(&std::fs::read_to_string(spec_path)?)?;
        let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);

        println!("Verifying {} against {}", path.display(), spec_path.display());

        // Compare each file
        let expected_types = generate_types(&spec, &algebra).to_string();
        let actual_types = std::fs::read_to_string(path.join("types.rs"))?;

        if normalize_whitespace(&expected_types) != normalize_whitespace(&actual_types) {
            println!("  types.rs: MISMATCH");
            return Err(anyhow!("Generated code doesn't match specification"));
        }

        println!("  types.rs: OK");

        // ... check other files
    }

    // Run the verification tests
    println!("Running verification tests...");
    let status = std::process::Command::new("cargo")
        .args(["test", "--", "--test-threads=1", "::verification::"])
        .status()?;

    if status.success() {
        println!("Verification PASSED");
        Ok(())
    } else {
        Err(anyhow!("Verification FAILED"))
    }
}
```

### Products Command

```rust
fn list_products(
    spec_path: &Path,
    type_filter: Option<&str>,
    product_filter: Option<&str>,
    table_format: bool,
) -> Result<()> {
    let spec = parse_spec(&std::fs::read_to_string(spec_path)?)?;
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);

    if table_format {
        println!("| A | B | Product | Result |");
        println!("|---|---|---------|--------|");
    }

    for a in &spec.types {
        if let Some(filter) = type_filter {
            if a.name != filter {
                continue;
            }
        }

        for b in &spec.types {
            // Geometric product
            if product_filter.is_none() || product_filter == Some("geometric") {
                if let Some(output) = find_geometric_output(&spec, a, b, &algebra) {
                    if table_format {
                        println!("| {} | {} | geometric | {} |", a.name, b.name, output.name);
                    } else {
                        println!("geometric({}, {}) -> {}", a.name, b.name, output.name);
                    }
                }
            }

            // Outer product
            if product_filter.is_none() || product_filter == Some("outer") {
                if let Some(output) = find_outer_output(&spec, a, b, &algebra) {
                    if table_format {
                        println!("| {} | {} | outer | {} |", a.name, b.name, output.name);
                    } else {
                        println!("outer({}, {}) -> {}", a.name, b.name, output.name);
                    }
                }
            }

            // ... other products
        }
    }

    Ok(())
}
```

## Build Integration

### Makefile Target

```makefile
# Makefile

ALGEBRAS := euclidean2 euclidean3 pga2 pga3 cga3

.PHONY: generate
generate:
	@for alg in $(ALGEBRAS); do \
		clifford-codegen generate algebras/$$alg.toml -o src/generated/$$alg/ --force; \
	done

.PHONY: verify-generated
verify-generated:
	@for alg in $(ALGEBRAS); do \
		clifford-codegen verify src/generated/$$alg/ --spec algebras/$$alg.toml || exit 1; \
	done
```

### Just Recipe

```just
# justfile

# Regenerate all algebras
generate:
    for spec in algebras/*.toml; do \
        clifford-codegen generate "$spec" --force
    done

# Verify generated code
verify:
    for spec in algebras/*.toml; do \
        name=$(basename "$spec" .toml)
        clifford-codegen verify "src/generated/$name" --spec "$spec"
    done
```

### CI Workflow

```yaml
# .github/workflows/verify-generated.yml

name: Verify Generated Code

on:
  push:
    paths:
      - 'algebras/*.toml'
      - 'src/generated/**'
  pull_request:
    paths:
      - 'algebras/*.toml'
      - 'src/generated/**'

jobs:
  verify:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install Rust
        uses: dtolnay/rust-action@stable

      - name: Install clifford-codegen
        run: cargo install --path crates/clifford-codegen

      - name: Verify generated code matches specs
        run: |
          for spec in algebras/*.toml; do
            name=$(basename "$spec" .toml)
            clifford-codegen verify "src/generated/$name" --spec "$spec"
          done

      - name: Regenerate and check for changes
        run: |
          clifford-codegen generate algebras/*.toml --force
          git diff --exit-code src/generated/
```

### Pre-commit Hook

```bash
#!/bin/bash
# .git/hooks/pre-commit

# Check if any algebra specs changed
if git diff --cached --name-only | grep -q "^algebras/.*\.toml$"; then
    echo "Algebra specifications changed. Verifying generated code..."

    # Regenerate and check for uncommitted changes
    clifford-codegen generate algebras/*.toml --force

    if ! git diff --quiet src/generated/; then
        echo "ERROR: Generated code out of sync with specifications."
        echo "Run 'make generate' and commit the changes."
        exit 1
    fi
fi
```

## Library API

For programmatic use:

```rust
use clifford_codegen::{parse_spec, generate_algebra, GenerateOptions};

fn main() -> Result<()> {
    let spec = parse_spec(r#"
        [algebra]
        name = "custom"

        [signature]
        positive = ["e1", "e2"]

        [types.Vector]
        grades = [1]
    "#)?;

    let options = GenerateOptions::default()
        .with_serde(true)
        .with_tests(true);

    let generated = generate_algebra(&spec, &options)?;

    // generated.types contains the types.rs TokenStream
    // generated.products contains the products.rs TokenStream
    // etc.

    Ok(())
}
```

## Installation

### From crates.io

```bash
cargo install clifford-codegen
```

### From source

```bash
git clone https://github.com/DevonMorris/clifford
cd clifford
cargo install --path crates/clifford-codegen
```

### As dev-dependency

```toml
# Cargo.toml
[dev-dependencies]
clifford-codegen = "0.1"
```

## Deliverables

- [ ] CLI argument parsing (clap)
- [ ] `generate` command implementation
- [ ] `verify` command implementation
- [ ] `products` command implementation
- [ ] `blades` command implementation
- [ ] rustfmt integration
- [ ] Library API for programmatic use
- [ ] Makefile/justfile examples
- [ ] CI workflow
- [ ] Pre-commit hook
- [ ] Installation documentation

## Dependencies

- `clap` - Command-line parsing
- `anyhow` - Error handling
- `proc-macro2` - Token manipulation
- `quote` - Quasi-quoting
- `toml` - Spec parsing

## Success Criteria

1. CLI generates correct, formatted code
2. Verification catches spec/code mismatches
3. CI prevents out-of-sync generated code
4. Library API enables custom tooling
5. Generated code compiles without modification
