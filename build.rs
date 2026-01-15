//! Build script that automatically regenerates all algebras from TOML specifications.
//!
//! This ensures generated files are always up-to-date with the codegen tool and TOML specs.
//! Never edit files in `src/specialized/*/generated/` - they will be overwritten on every build.

use std::path::Path;

use clifford_codegen::algebra::{Algebra, ProductTable};
use clifford_codegen::codegen::{
    ConversionsGenerator, TraitsGenerator, TypeGenerator, format_tokens,
};
use clifford_codegen::spec::parse_spec;

/// Algebra configurations: (name, toml_path, output_dir)
const ALGEBRAS: &[(&str, &str, &str)] = &[
    (
        "euclidean2",
        "algebras/euclidean2.toml",
        "src/specialized/euclidean/dim2/generated",
    ),
    (
        "euclidean3",
        "algebras/euclidean3.toml",
        "src/specialized/euclidean/dim3/generated",
    ),
    (
        "projective2",
        "algebras/projective2.toml",
        "src/specialized/projective/dim2/generated",
    ),
    (
        "projective3",
        "algebras/projective3.toml",
        "src/specialized/projective/dim3/generated",
    ),
    (
        "hyperbolic",
        "algebras/hyperbolic.toml",
        "src/specialized/hyperbolic/generated",
    ),
    (
        "complex",
        "algebras/complex.toml",
        "src/specialized/complex/generated",
    ),
    (
        "dual",
        "algebras/dual.toml",
        "src/specialized/dual/generated",
    ),
    (
        "quaternion",
        "algebras/quaternion.toml",
        "src/specialized/quaternion/generated",
    ),
    (
        "minkowski2",
        "algebras/minkowski2.toml",
        "src/specialized/minkowski/dim2/generated",
    ),
    (
        "dualquat",
        "algebras/dualquat.toml",
        "src/specialized/dualquat/generated",
    ),
    (
        "elliptic2",
        "algebras/elliptic2.toml",
        "src/specialized/elliptic/dim2/generated",
    ),
];

fn main() {
    for (name, toml_path, output_dir) in ALGEBRAS {
        generate_algebra(name, toml_path, output_dir);
    }
}

/// Generates code for a single algebra.
fn generate_algebra(name: &str, toml_path: &str, output_dir: &str) {
    eprintln!("Regenerating {}...", name);

    // Parse specification
    let spec_content = std::fs::read_to_string(toml_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {}", toml_path, e));
    let spec = parse_spec(&spec_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {}", toml_path, e));

    // Build algebra and product table
    let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
    let table = ProductTable::new(&algebra);

    // Create generators
    let type_gen = TypeGenerator::new(&spec, &algebra);
    let traits_gen = TraitsGenerator::new(&spec, &algebra, table);
    let conversions_gen = ConversionsGenerator::new(&spec, &algebra);

    // Create output directory
    let output_path = Path::new(output_dir);
    std::fs::create_dir_all(output_path)
        .unwrap_or_else(|e| panic!("Failed to create {}: {}", output_dir, e));

    // Generate and write mod.rs
    let mod_content = generate_mod_file(&spec.name);
    write_file(&output_path.join("mod.rs"), &format_tokens(&mod_content));

    // Generate and write types.rs
    let types_content = type_gen.generate_types_file();
    write_file(
        &output_path.join("types.rs"),
        &format_tokens(&types_content),
    );

    // Generate and write conversions.rs
    let conversions_content = conversions_gen.generate_conversions_file();
    write_file(
        &output_path.join("conversions.rs"),
        &format_tokens(&conversions_content),
    );

    // Generate and write traits.rs (special: returns tokens + formatted tests string)
    let (traits_tokens, traits_tests) = traits_gen.generate_traits_file();
    let traits_content = format!("{}{}", format_tokens(&traits_tokens), traits_tests);
    write_file(&output_path.join("traits.rs"), &traits_content);
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
    std::fs::write(path, content)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", path.display(), e));

    // Format the generated file with rustfmt
    let status = std::process::Command::new("rustfmt")
        .arg("--edition")
        .arg("2024")
        .arg(path)
        .status();

    match status {
        Ok(s) if s.success() => {}
        Ok(s) => eprintln!(
            "rustfmt failed for {}: exit code {:?}",
            path.display(),
            s.code()
        ),
        Err(e) => eprintln!("Failed to run rustfmt for {}: {}", path.display(), e),
    }
}
