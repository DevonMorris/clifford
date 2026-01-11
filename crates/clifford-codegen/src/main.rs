//! CLI tool for generating optimized geometric algebra code.
//!
//! This tool generates Rust code from algebra specification files.
//!
//! # Usage
//!
//! ```bash
//! # Generate code from a specification
//! clifford-codegen generate algebras/euclidean3.toml -o src/generated/
//!
//! # Verify generated code matches specification
//! clifford-codegen verify src/generated/euclidean3/ --spec algebras/euclidean3.toml
//!
//! # List products for an algebra
//! clifford-codegen products algebras/pga3.toml --type Motor
//!
//! # List blades
//! clifford-codegen blades algebras/euclidean3.toml
//! ```

use std::path::PathBuf;

use anyhow::{Result, anyhow};
use clap::{Parser, Subcommand};

use clifford_codegen::algebra::{Algebra, ProductTable};
use clifford_codegen::codegen::{
    ConstraintGenerator, ConversionsGenerator, ProductGenerator, TraitsGenerator, TypeGenerator,
    format_tokens,
};
use clifford_codegen::discovery::{discover_entities, generate_toml_template};
use clifford_codegen::spec::parse_spec;

/// Generate optimized geometric algebra code from specifications.
#[derive(Parser)]
#[command(name = "clifford-codegen")]
#[command(version, about)]
struct Cli {
    /// Subcommand to execute.
    #[command(subcommand)]
    command: Commands,
}

/// Available subcommands.
#[derive(Subcommand)]
enum Commands {
    /// Generate Rust code from algebra specifications.
    Generate {
        /// Specification file(s).
        #[arg(required = true)]
        specs: Vec<PathBuf>,

        /// Output directory.
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Include nalgebra conversions.
        #[arg(long)]
        nalgebra: bool,

        /// Include serde derives.
        #[arg(long)]
        serde: bool,

        /// Skip test generation.
        #[arg(long)]
        no_tests: bool,

        /// Skip Arbitrary implementations.
        #[arg(long)]
        no_arbitrary: bool,

        /// Show what would be generated.
        #[arg(long)]
        dry_run: bool,

        /// Overwrite existing files.
        #[arg(long)]
        force: bool,

        /// Verbose output.
        #[arg(short, long)]
        verbose: bool,
    },

    /// Verify generated code matches specification.
    Verify {
        /// Path to generated module.
        path: PathBuf,

        /// Original specification file.
        #[arg(long)]
        spec: Option<PathBuf>,
    },

    /// List products for an algebra.
    Products {
        /// Specification file.
        spec: PathBuf,

        /// Filter by type.
        #[arg(long, name = "type")]
        type_name: Option<String>,

        /// Filter by product kind (geometric, outer, inner, sandwich).
        #[arg(long)]
        product: Option<String>,

        /// Output as markdown table.
        #[arg(long)]
        table: bool,
    },

    /// List blades and their indices.
    Blades {
        /// Specification file.
        spec: PathBuf,

        /// Filter by grade.
        #[arg(long)]
        grade: Option<usize>,

        /// Output format: table, json.
        #[arg(long, default_value = "table")]
        format: String,
    },

    /// Discover valid entities for an algebra signature.
    ///
    /// Analyzes which grade combinations satisfy geometric constraints
    /// and generates a TOML template for further customization.
    Discover {
        /// Algebra signature as "p,q,r" (positive, negative, zero dimensions).
        /// Examples: "3,0,0" for Euclidean 3D, "3,0,1" for PGA 3D, "4,1,0" for CGA 3D.
        signature: String,

        /// Output file for the TOML template. If not specified, prints to stdout.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

/// Options for code generation.
#[allow(dead_code)]
struct GenerateOptions {
    /// Output directory.
    output: Option<PathBuf>,
    /// Include nalgebra conversions.
    nalgebra: bool,
    /// Include serde derives.
    serde: bool,
    /// Generate tests.
    tests: bool,
    /// Generate Arbitrary implementations.
    arbitrary: bool,
    /// Dry run mode.
    dry_run: bool,
    /// Force overwrite.
    force: bool,
    /// Verbose output.
    verbose: bool,
}

fn main() -> Result<()> {
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
            let options = GenerateOptions {
                output,
                nalgebra,
                serde,
                tests: !no_tests,
                arbitrary: !no_arbitrary,
                dry_run,
                force,
                verbose,
            };

            for spec_path in specs {
                generate(&spec_path, &options)?;
            }
        }
        Commands::Verify { path, spec } => {
            verify(&path, spec.as_deref())?;
        }
        Commands::Products {
            spec,
            type_name,
            product,
            table,
        } => {
            list_products(&spec, type_name.as_deref(), product.as_deref(), table)?;
        }
        Commands::Blades {
            spec,
            grade,
            format,
        } => {
            list_blades(&spec, grade, &format)?;
        }
        Commands::Discover { signature, output } => {
            discover_algebra(&signature, output.as_deref())?;
        }
    }

    Ok(())
}

/// Generates code for a single algebra specification.
fn generate(spec_path: &std::path::Path, options: &GenerateOptions) -> Result<()> {
    // Parse specification
    let spec_content = std::fs::read_to_string(spec_path)?;
    let spec = parse_spec(&spec_content)?;

    if options.verbose {
        println!("Generating algebra: {}", spec.name);
        println!(
            "  Signature: ({}, {}, {})",
            spec.signature.p, spec.signature.q, spec.signature.r
        );
        println!("  Types: {}", spec.types.len());
    }

    // Determine output directory
    let output_dir = options
        .output
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("src/generated/{}", spec.name)));

    // Check if output exists
    if output_dir.exists() && !options.force && !options.dry_run {
        return Err(anyhow!(
            "Output directory exists. Use --force to overwrite."
        ));
    }

    // Build algebra
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
    let table = ProductTable::new(&algebra);

    // Create generators
    let type_gen = TypeGenerator::new(&spec, &algebra);
    let product_gen = ProductGenerator::new(&spec, &algebra, table.clone());
    let constraint_gen = ConstraintGenerator::new(&spec, &algebra);
    let traits_gen = TraitsGenerator::new(&spec, &algebra, table);
    let conversions_gen = ConversionsGenerator::new(&spec, &algebra);

    // Generate all files
    let files = vec![
        ("types.rs", type_gen.generate_types_file()),
        ("products.rs", product_gen.generate_products_file()),
        ("constrained.rs", constraint_gen.generate_constraints_file()),
        ("traits.rs", traits_gen.generate_traits_file()),
        (
            "conversions.rs",
            conversions_gen.generate_conversions_file(),
        ),
    ];

    // Generate mod.rs
    let mod_content = generate_mod_file(&spec, options);

    if options.dry_run {
        println!("Would generate in {}:", output_dir.display());
        println!("  mod.rs");
        for (name, _) in &files {
            println!("  {}", name);
        }
        return Ok(());
    }

    // Create output directory
    std::fs::create_dir_all(&output_dir)?;

    // Write mod.rs
    let mod_path = output_dir.join("mod.rs");
    let formatted_mod = format_tokens(&mod_content);
    std::fs::write(&mod_path, formatted_mod)?;
    if options.verbose {
        println!("  Wrote: {}", mod_path.display());
    }

    // Write other files
    for (name, tokens) in files {
        let path = output_dir.join(name);
        let formatted = format_tokens(&tokens);
        std::fs::write(&path, formatted)?;
        if options.verbose {
            println!("  Wrote: {}", path.display());
        }
    }

    println!("Generated {} in {}", spec.name, output_dir.display());

    Ok(())
}

/// Generates the mod.rs file content.
fn generate_mod_file(
    spec: &clifford_codegen::spec::AlgebraSpec,
    options: &GenerateOptions,
) -> proc_macro2::TokenStream {
    use quote::quote;

    let name = &spec.name;
    let doc = format!("Generated code for the {} algebra.", name);

    let mut mods = vec![
        quote! { pub mod types; },
        quote! { pub mod products; },
        quote! { pub mod constrained; },
        quote! { pub mod traits; },
        quote! { pub mod conversions; },
    ];

    if options.nalgebra {
        mods.push(quote! { pub mod nalgebra; });
    }

    quote! {
        #![doc = #doc]

        #(#mods)*
    }
}

/// Verifies that generated code matches the specification.
fn verify(path: &std::path::Path, spec_path: Option<&std::path::Path>) -> Result<()> {
    // Check that the path exists
    let mod_rs = path.join("mod.rs");
    if !mod_rs.exists() {
        return Err(anyhow!("Not a generated module: {}", path.display()));
    }

    println!("Verifying module at: {}", path.display());

    // If spec provided, regenerate and compare
    if let Some(spec_path) = spec_path {
        let spec_content = std::fs::read_to_string(spec_path)?;
        let spec = parse_spec(&spec_content)?;
        let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
        let table = ProductTable::new(&algebra);

        println!("  Comparing against: {}", spec_path.display());

        // Check types.rs
        let type_gen = TypeGenerator::new(&spec, &algebra);
        let expected_types = format_tokens(&type_gen.generate_types_file());
        let types_path = path.join("types.rs");

        if types_path.exists() {
            let actual_types = std::fs::read_to_string(&types_path)?;
            if normalize_whitespace(&expected_types) != normalize_whitespace(&actual_types) {
                println!("  types.rs: MISMATCH");
                return Err(anyhow!("Generated code doesn't match specification"));
            }
            println!("  types.rs: OK");
        }

        // Check products.rs
        let product_gen = ProductGenerator::new(&spec, &algebra, table);
        let expected_products = format_tokens(&product_gen.generate_products_file());
        let products_path = path.join("products.rs");

        if products_path.exists() {
            let actual_products = std::fs::read_to_string(&products_path)?;
            if normalize_whitespace(&expected_products) != normalize_whitespace(&actual_products) {
                println!("  products.rs: MISMATCH");
                return Err(anyhow!("Generated code doesn't match specification"));
            }
            println!("  products.rs: OK");
        }
    }

    println!("Verification complete.");
    Ok(())
}

/// Normalizes whitespace for comparison.
fn normalize_whitespace(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Lists products for an algebra.
fn list_products(
    spec_path: &std::path::Path,
    type_filter: Option<&str>,
    product_filter: Option<&str>,
    table_format: bool,
) -> Result<()> {
    let spec_content = std::fs::read_to_string(spec_path)?;
    let spec = parse_spec(&spec_content)?;
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);

    if table_format {
        println!("| A | B | Product | Result Type |");
        println!("|---|---|---------|-------------|");
    }

    // Collect type names
    let type_names: Vec<_> = spec.types.iter().map(|t| t.name.as_str()).collect();

    for type_a in &spec.types {
        if let Some(filter) = type_filter {
            if type_a.name != filter {
                continue;
            }
        }

        for type_b in &spec.types {
            // Check geometric product
            if product_filter.is_none() || product_filter == Some("geometric") {
                if let Some(result) =
                    find_product_result(&type_a.grades, &type_b.grades, &algebra, &type_names)
                {
                    if table_format {
                        println!(
                            "| {} | {} | geometric | {} |",
                            type_a.name, type_b.name, result
                        );
                    } else {
                        println!("geometric({}, {}) -> {}", type_a.name, type_b.name, result);
                    }
                }
            }

            // Check outer product
            if product_filter.is_none() || product_filter == Some("outer") {
                if let Some(result) =
                    find_outer_result(&type_a.grades, &type_b.grades, &algebra, &type_names)
                {
                    if table_format {
                        println!("| {} | {} | outer | {} |", type_a.name, type_b.name, result);
                    } else {
                        println!("outer({}, {}) -> {}", type_a.name, type_b.name, result);
                    }
                }
            }
        }
    }

    Ok(())
}

/// Finds the result type for a geometric product.
fn find_product_result(
    grades_a: &[usize],
    grades_b: &[usize],
    _algebra: &Algebra,
    _type_names: &[&str],
) -> Option<String> {
    // For geometric product, result grades are |grade_a - grade_b| to grade_a + grade_b
    let mut result_grades = std::collections::BTreeSet::new();
    for &ga in grades_a {
        for &gb in grades_b {
            let min_grade = ga.abs_diff(gb);
            let max_grade = ga + gb;
            for g in min_grade..=max_grade {
                if (g + min_grade) % 2 == 0 {
                    // Only grades with same parity difference
                    result_grades.insert(g);
                }
            }
        }
    }

    if result_grades.is_empty() {
        None
    } else {
        Some(format!(
            "Grade({:?})",
            result_grades.iter().collect::<Vec<_>>()
        ))
    }
}

/// Finds the result type for an outer product.
fn find_outer_result(
    grades_a: &[usize],
    grades_b: &[usize],
    _algebra: &Algebra,
    _type_names: &[&str],
) -> Option<String> {
    // For outer product, result grade is grade_a + grade_b
    let mut result_grades = std::collections::BTreeSet::new();
    for &ga in grades_a {
        for &gb in grades_b {
            result_grades.insert(ga + gb);
        }
    }

    if result_grades.is_empty() {
        None
    } else {
        Some(format!(
            "Grade({:?})",
            result_grades.iter().collect::<Vec<_>>()
        ))
    }
}

/// Lists blades for an algebra.
fn list_blades(
    spec_path: &std::path::Path,
    grade_filter: Option<usize>,
    format: &str,
) -> Result<()> {
    let spec_content = std::fs::read_to_string(spec_path)?;
    let spec = parse_spec(&spec_content)?;
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);

    let dim = algebra.dim();
    let num_blades = 1usize << dim;

    if format == "json" {
        println!("[");
    } else {
        println!("| Index | Grade | Blade | Metric |");
        println!("|-------|-------|-------|--------|");
    }

    let mut first = true;
    for index in 0..num_blades {
        let grade = (index as u32).count_ones() as usize;

        if let Some(filter) = grade_filter {
            if grade != filter {
                continue;
            }
        }

        let blade_name_str = compute_blade_name(index, &spec);
        let metric = blade_metric(index, &algebra);

        if format == "json" {
            if !first {
                println!(",");
            }
            first = false;
            print!(
                "  {{ \"index\": {}, \"grade\": {}, \"name\": \"{}\", \"metric\": {} }}",
                index, grade, blade_name_str, metric
            );
        } else {
            println!(
                "| {} | {} | {} | {} |",
                index, grade, blade_name_str, metric
            );
        }
    }

    if format == "json" {
        println!("\n]");
    }

    Ok(())
}

/// Computes the blade name from its index.
fn compute_blade_name(index: usize, spec: &clifford_codegen::spec::AlgebraSpec) -> String {
    if index == 0 {
        return "1".to_string();
    }

    let mut parts = Vec::new();
    for (i, basis) in spec.signature.basis.iter().enumerate() {
        if (index >> i) & 1 == 1 {
            parts.push(basis.name.as_str());
        }
    }

    parts.join("")
}

/// Computes the metric of a blade.
fn blade_metric(index: usize, algebra: &Algebra) -> i32 {
    let mut result = 1i32;
    let dim = algebra.dim();

    for i in 0..dim {
        if (index >> i) & 1 == 1 {
            result *= i32::from(algebra.metric(i));
        }
    }

    result
}

/// Discovers valid entities for an algebra and generates a TOML template.
fn discover_algebra(signature_str: &str, output_path: Option<&std::path::Path>) -> Result<()> {
    // Parse signature "p,q,r"
    let parts: Vec<&str> = signature_str.split(',').collect();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(anyhow!(
            "Invalid signature format. Expected 'p,q' or 'p,q,r', got '{}'",
            signature_str
        ));
    }

    let p: usize = parts[0]
        .trim()
        .parse()
        .map_err(|_| anyhow!("Invalid positive dimension: '{}'", parts[0]))?;
    let q: usize = parts[1]
        .trim()
        .parse()
        .map_err(|_| anyhow!("Invalid negative dimension: '{}'", parts[1]))?;
    let r: usize = if parts.len() > 2 {
        parts[2]
            .trim()
            .parse()
            .map_err(|_| anyhow!("Invalid zero dimension: '{}'", parts[2]))?
    } else {
        0
    };

    // Create algebra
    let algebra = Algebra::new(p, q, r);
    let (p, q, r) = algebra.signature();

    eprintln!(
        "Discovering entities for Cl({},{},{})...",
        p, q, r
    );

    // Discover entities
    let entities = discover_entities(&algebra);

    eprintln!(
        "Found {} valid grade combinations that satisfy geometric constraints.",
        entities.len()
    );

    // Generate TOML template
    if let Some(path) = output_path {
        let mut file = std::fs::File::create(path)?;
        generate_toml_template(&algebra, &entities, &mut file)?;
        eprintln!("Wrote template to: {}", path.display());
    } else {
        let mut stdout = std::io::stdout();
        generate_toml_template(&algebra, &entities, &mut stdout)?;
    }

    Ok(())
}
