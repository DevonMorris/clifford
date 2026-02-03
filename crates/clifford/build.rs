//! Build script for clifford.
//!
//! When the `codegen` feature is enabled, this regenerates all algebras from TOML specs.
//! When `codegen` is disabled (default for crates.io users), it verifies generated files exist.

fn main() {
    // Check if codegen feature is enabled
    #[cfg(feature = "codegen")]
    {
        regenerate_algebras();
    }

    #[cfg(not(feature = "codegen"))]
    {
        verify_generated_files_exist();
    }
}

/// Verify that all generated files exist (for crates.io users without codegen).
#[cfg(not(feature = "codegen"))]
fn verify_generated_files_exist() {
    use std::path::Path;

    const GENERATED_DIRS: &[&str] = &[
        "src/specialized/euclidean/dim2/generated",
        "src/specialized/euclidean/dim3/generated",
        "src/specialized/projective/dim2/generated",
        "src/specialized/projective/dim3/generated",
        "src/specialized/complex/generated",
        "src/specialized/quaternion/generated",
        "src/specialized/dual/generated",
        "src/specialized/hyperbolic/generated",
        "src/specialized/minkowski/dim2/generated",
        "src/specialized/dualquat/generated",
        "src/specialized/elliptic/dim2/generated",
        "src/specialized/hyperbolic/dim2/generated",
        "src/specialized/minkowski/dim3/generated",
        "src/specialized/conformal/dim3/generated",
        "src/specialized/conformal/dim2/generated",
    ];

    for dir in GENERATED_DIRS {
        let types_path = Path::new(dir).join("types.rs");
        if !types_path.exists() {
            panic!(
                "Generated file {} does not exist.\n\
                 \n\
                 If you're building from git, enable the `codegen` feature:\n\
                 \n\
                 cargo build --features codegen\n\
                 \n\
                 This requires a Symbolica license for development.\n\
                 \n\
                 If you're building from crates.io, this is a packaging bug - please report it.",
                types_path.display()
            );
        }
    }
}

/// Regenerate all algebras from TOML specs (requires codegen feature).
#[cfg(feature = "codegen")]
fn regenerate_algebras() {
    use std::path::Path;
    use std::time::Instant;

    use clifford_codegen::algebra::{Algebra, ProductTable};
    use clifford_codegen::codegen::{
        ConversionsGenerator, TraitsGenerator, TypeGenerator, format_tokens,
    };
    use clifford_codegen::spec::parse_spec;
    use proc_macro2::TokenStream;
    use quote::quote;

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
        use_groebner: bool,
    }

    impl Default for AlgebraOptions {
        fn default() -> Self {
            Self { use_groebner: true }
        }
    }

    /// Configuration for a single algebra to generate.
    struct AlgebraConfig {
        name: &'static str,
        toml_path: &'static str,
        output_dir: &'static str,
        options: AlgebraOptions,
    }

    impl AlgebraConfig {
        const fn new(
            name: &'static str,
            toml_path: &'static str,
            output_dir: &'static str,
        ) -> Self {
            Self {
                name,
                toml_path,
                output_dir,
                options: AlgebraOptions { use_groebner: true },
            }
        }

        const fn with_no_groebner(mut self) -> Self {
            self.options.use_groebner = false;
            self
        }
    }

    const ALGEBRAS: &[AlgebraConfig] = &[
        AlgebraConfig::new(
            "euclidean2",
            "../../algebras/euclidean2.toml",
            "src/specialized/euclidean/dim2/generated",
        ),
        AlgebraConfig::new(
            "euclidean3",
            "../../algebras/euclidean3.toml",
            "src/specialized/euclidean/dim3/generated",
        ),
        AlgebraConfig::new(
            "projective2",
            "../../algebras/projective2.toml",
            "src/specialized/projective/dim2/generated",
        ),
        AlgebraConfig::new(
            "projective3",
            "../../algebras/projective3.toml",
            "src/specialized/projective/dim3/generated",
        ),
        AlgebraConfig::new(
            "complex",
            "../../algebras/complex.toml",
            "src/specialized/complex/generated",
        ),
        AlgebraConfig::new(
            "quaternion",
            "../../algebras/quaternion.toml",
            "src/specialized/quaternion/generated",
        ),
        AlgebraConfig::new(
            "dual",
            "../../algebras/dual.toml",
            "src/specialized/dual/generated",
        ),
        AlgebraConfig::new(
            "hyperbolic",
            "../../algebras/hyperbolic.toml",
            "src/specialized/hyperbolic/generated",
        ),
        AlgebraConfig::new(
            "minkowski2",
            "../../algebras/minkowski2.toml",
            "src/specialized/minkowski/dim2/generated",
        ),
        AlgebraConfig::new(
            "dualquat",
            "../../algebras/dualquat.toml",
            "src/specialized/dualquat/generated",
        ),
        AlgebraConfig::new(
            "elliptic2",
            "../../algebras/elliptic2.toml",
            "src/specialized/elliptic/dim2/generated",
        ),
        AlgebraConfig::new(
            "hyperbolic2",
            "../../algebras/hyperbolic2.toml",
            "src/specialized/hyperbolic/dim2/generated",
        ),
        AlgebraConfig::new(
            "minkowski3",
            "../../algebras/minkowski3.toml",
            "src/specialized/minkowski/dim3/generated",
        ),
        AlgebraConfig::new(
            "conformal3",
            "../../algebras/conformal3.toml",
            "src/specialized/conformal/dim3/generated",
        )
        .with_no_groebner(),
        AlgebraConfig::new(
            "conformal2",
            "../../algebras/conformal2.toml",
            "src/specialized/conformal/dim2/generated",
        )
        .with_no_groebner(),
    ];

    fn generate_mod_file(name: &str) -> TokenStream {
        let doc = format!("Generated code for the {} algebra.", name);
        quote! {
            #![doc = #doc]

            pub mod types;
            pub mod traits;
            pub mod conversions;
        }
    }

    fn write_file(path: &Path, content: &str) {
        use std::fs;
        use std::io::Write;

        log_debug!("    Writing {} ({} bytes)", path.display(), content.len());
        let mut file = fs::File::create(path)
            .unwrap_or_else(|e| panic!("Failed to create {}: {}", path.display(), e));
        file.write_all(content.as_bytes())
            .unwrap_or_else(|e| panic!("Failed to write {}: {}", path.display(), e));
    }

    fn format_file(path: &Path) {
        use std::process::Command;

        log_debug!("    Formatted {}", path.display());
        let status = Command::new("rustfmt")
            .arg("--edition")
            .arg("2024")
            .arg(path)
            .status();

        if let Err(e) = status {
            eprintln!("Warning: rustfmt failed for {}: {}", path.display(), e);
        }
    }

    fn generate_algebra(config: &AlgebraConfig) {
        let start = Instant::now();

        log_debug!("  Parsing spec from {}", config.toml_path);
        let spec_content = std::fs::read_to_string(config.toml_path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {}", config.toml_path, e));
        let spec = parse_spec(&spec_content)
            .unwrap_or_else(|e| panic!("Failed to parse {}: {}", config.toml_path, e));

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

        let type_gen = TypeGenerator::new(&spec, &algebra);
        let traits_gen =
            TraitsGenerator::with_options(&spec, &algebra, table, config.options.use_groebner);
        let conversions_gen = ConversionsGenerator::new(&spec, &algebra);

        let output_path = Path::new(config.output_dir);
        std::fs::create_dir_all(output_path)
            .unwrap_or_else(|e| panic!("Failed to create {}: {}", config.output_dir, e));

        let mod_path = output_path.join("mod.rs");
        let types_path = output_path.join("types.rs");
        let conversions_path = output_path.join("conversions.rs");
        let traits_path = output_path.join("traits.rs");

        log_debug!("  Generating mod.rs");
        let mod_content = generate_mod_file(&spec.name);
        write_file(&mod_path, &format_tokens(&mod_content));

        log_debug!("  Generating types.rs");
        let types_content = type_gen.generate_types_file();
        write_file(&types_path, &format_tokens(&types_content));

        log_debug!("  Generating conversions.rs");
        let conversions_content = conversions_gen.generate_conversions_file();
        write_file(&conversions_path, &format_tokens(&conversions_content));

        log_debug!("  Generating traits.rs");
        let (traits_content, test_content) = traits_gen.generate_traits_file();
        let full_traits = format!("{}{}", format_tokens(&traits_content), test_content);
        write_file(&traits_path, &full_traits);

        format_file(&mod_path);
        format_file(&types_path);
        format_file(&conversions_path);
        format_file(&traits_path);

        let elapsed = start.elapsed();
        log_debug!(
            "  Completed {} in {:.2}ms",
            config.name,
            elapsed.as_secs_f64() * 1000.0
        );
    }

    // Main codegen logic
    let total_start = Instant::now();
    log_info!("Starting code generation for {} algebras", ALGEBRAS.len());

    for config in ALGEBRAS {
        println!("cargo::rerun-if-changed={}", config.toml_path);
    }
    println!("cargo::rerun-if-changed=../clifford-codegen/src");

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
