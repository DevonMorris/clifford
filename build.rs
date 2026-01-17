//! Build script that automatically regenerates all algebras from TOML specifications.
//!
//! This ensures generated files are always up-to-date with the codegen tool and TOML specs.
//! Never edit files in `src/specialized/*/generated/` - they will be overwritten on every build.

use std::path::Path;
use std::time::Instant;

use clifford_codegen::algebra::{Algebra, ProductTable};
use clifford_codegen::codegen::{
    ConversionsGenerator, TraitsGenerator, TypeGenerator, format_tokens,
};
use clifford_codegen::spec::parse_spec;

/// Logs a message that will be visible during `cargo build`.
macro_rules! log_info {
    ($($arg:tt)*) => {
        println!("cargo::warning=[clifford] {}", format!($($arg)*))
    };
}

/// Logs a debug message (only visible with `cargo build -vv`).
macro_rules! log_debug {
    ($($arg:tt)*) => {
        eprintln!("[clifford build] {}", format!($($arg)*))
    };
}

/// Options for algebra code generation.
struct AlgebraOptions {
    /// Whether to use Groebner basis simplification for constrained types.
    ///
    /// When true (default), the generator uses Groebner basis reduction to simplify
    /// expressions for constrained types. This can be disabled for algebras where
    /// Groebner computation is too expensive (e.g., high-dimensional CGA).
    use_groebner: bool,
}

impl Default for AlgebraOptions {
    fn default() -> Self {
        Self { use_groebner: true }
    }
}

/// Configuration for a single algebra to generate.
struct AlgebraConfig {
    /// Display name for the algebra (used in logging).
    name: &'static str,
    /// Path to the TOML specification file.
    toml_path: &'static str,
    /// Output directory for generated code.
    output_dir: &'static str,
    /// Code generation options.
    options: AlgebraOptions,
}

impl AlgebraConfig {
    /// Creates a new algebra configuration with default options.
    const fn new(name: &'static str, toml_path: &'static str, output_dir: &'static str) -> Self {
        Self {
            name,
            toml_path,
            output_dir,
            options: AlgebraOptions { use_groebner: true },
        }
    }

    /// Disables Groebner basis simplification for this algebra.
    const fn with_no_groebner(mut self) -> Self {
        self.options.use_groebner = false;
        self
    }
}

/// Algebra configurations.
const ALGEBRAS: &[AlgebraConfig] = &[
    AlgebraConfig::new(
        "euclidean2",
        "algebras/euclidean2.toml",
        "src/specialized/euclidean/dim2/generated",
    ),
    AlgebraConfig::new(
        "euclidean3",
        "algebras/euclidean3.toml",
        "src/specialized/euclidean/dim3/generated",
    ),
    AlgebraConfig::new(
        "projective2",
        "algebras/projective2.toml",
        "src/specialized/projective/dim2/generated",
    ),
    AlgebraConfig::new(
        "projective3",
        "algebras/projective3.toml",
        "src/specialized/projective/dim3/generated",
    ),
    AlgebraConfig::new(
        "hyperbolic",
        "algebras/hyperbolic.toml",
        "src/specialized/hyperbolic/generated",
    ),
    AlgebraConfig::new(
        "complex",
        "algebras/complex.toml",
        "src/specialized/complex/generated",
    ),
    AlgebraConfig::new(
        "dual",
        "algebras/dual.toml",
        "src/specialized/dual/generated",
    ),
    AlgebraConfig::new(
        "quaternion",
        "algebras/quaternion.toml",
        "src/specialized/quaternion/generated",
    ),
    AlgebraConfig::new(
        "minkowski2",
        "algebras/minkowski2.toml",
        "src/specialized/minkowski/dim2/generated",
    ),
    AlgebraConfig::new(
        "dualquat",
        "algebras/dualquat.toml",
        "src/specialized/dualquat/generated",
    ),
    AlgebraConfig::new(
        "elliptic2",
        "algebras/elliptic2.toml",
        "src/specialized/elliptic/dim2/generated",
    ),
    AlgebraConfig::new(
        "hyperbolic2",
        "algebras/hyperbolic2.toml",
        "src/specialized/hyperbolic/dim2/generated",
    ),
    AlgebraConfig::new(
        "minkowski3",
        "algebras/minkowski3.toml",
        "src/specialized/minkowski/dim3/generated",
    ),
    AlgebraConfig::new(
        "conformal3",
        "algebras/conformal3.toml",
        "src/specialized/conformal/dim3/generated",
    )
    .with_no_groebner(),
];

fn main() {
    let total_start = Instant::now();
    log_info!("Starting code generation for {} algebras", ALGEBRAS.len());

    // Register rerun-if-changed for all TOML files and the codegen crate
    for config in ALGEBRAS {
        println!("cargo::rerun-if-changed={}", config.toml_path);
    }
    println!("cargo::rerun-if-changed=crates/clifford-codegen/src");

    for (i, config) in ALGEBRAS.iter().enumerate() {
        log_debug!(
            "[{}/{}] Generating {}...",
            i + 1,
            ALGEBRAS.len(),
            config.name
        );
        generate_algebra(config);
    }

    let total_elapsed = total_start.elapsed();
    log_info!(
        "Code generation complete: {} algebras in {:.2}s",
        ALGEBRAS.len(),
        total_elapsed.as_secs_f64()
    );
}

/// Generates code for a single algebra.
fn generate_algebra(config: &AlgebraConfig) {
    let start = Instant::now();

    // Parse specification
    log_debug!("  Parsing spec from {}", config.toml_path);
    let spec_content = std::fs::read_to_string(config.toml_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {}", config.toml_path, e));
    let spec = parse_spec(&spec_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {}", config.toml_path, e));

    // Build algebra and product table
    log_debug!(
        "  Building algebra: signature ({}, {}, {}), {} types, groebner={}",
        spec.signature.p,
        spec.signature.q,
        spec.signature.r,
        spec.types.len(),
        config.options.use_groebner
    );
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
    let table = ProductTable::new(&algebra);

    // Create generators
    let type_gen = TypeGenerator::new(&spec, &algebra);
    let traits_gen =
        TraitsGenerator::with_options(&spec, &algebra, table, config.options.use_groebner);
    let conversions_gen = ConversionsGenerator::new(&spec, &algebra);

    // Create output directory
    let output_path = Path::new(config.output_dir);
    std::fs::create_dir_all(output_path)
        .unwrap_or_else(|e| panic!("Failed to create {}: {}", config.output_dir, e));

    // Generate and write mod.rs
    log_debug!("  Generating mod.rs");
    let mod_content = generate_mod_file(&spec.name);
    write_file(&output_path.join("mod.rs"), &format_tokens(&mod_content));

    // Generate and write types.rs
    log_debug!("  Generating types.rs");
    let types_content = type_gen.generate_types_file();
    write_file(
        &output_path.join("types.rs"),
        &format_tokens(&types_content),
    );

    // Generate and write conversions.rs
    log_debug!("  Generating conversions.rs");
    let conversions_content = conversions_gen.generate_conversions_file();
    write_file(
        &output_path.join("conversions.rs"),
        &format_tokens(&conversions_content),
    );

    // Generate and write traits.rs (special: returns tokens + formatted tests string)
    log_debug!("  Generating traits.rs");
    let (traits_tokens, traits_tests) = traits_gen.generate_traits_file();
    let traits_content = format!("{}{}", format_tokens(&traits_tokens), traits_tests);
    write_file(&output_path.join("traits.rs"), &traits_content);

    let elapsed = start.elapsed();
    log_debug!(
        "  Completed {} in {:.2}ms",
        config.name,
        elapsed.as_secs_f64() * 1000.0
    );
}

/// Generates the mod.rs file content.
fn generate_mod_file(name: &str) -> proc_macro2::TokenStream {
    use quote::quote;

    let doc = format!("Generated code for the {} algebra.", name);

    quote! {
        #![doc = #doc]

        pub mod types;
        pub mod traits;
        pub mod conversions;
    }
}

/// Writes content to a file and formats it with rustfmt.
fn write_file(path: &Path, content: &str) {
    log_debug!("    Writing {} ({} bytes)", path.display(), content.len());

    std::fs::write(path, content)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", path.display(), e));

    // Format the generated file with rustfmt
    let status = std::process::Command::new("rustfmt")
        .arg("--edition")
        .arg("2024")
        .arg(path)
        .status();

    match status {
        Ok(s) if s.success() => log_debug!("    Formatted {}", path.display()),
        Ok(s) => log_info!(
            "rustfmt failed for {}: exit code {:?}",
            path.display(),
            s.code()
        ),
        Err(e) => log_info!("Failed to run rustfmt for {}: {}", path.display(), e),
    }
}
