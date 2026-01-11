# Release Agent

You are handling release and publishing tasks for Clifford, a Rust geometric algebra library.

## Responsibilities

- Version management (semver)
- CHANGELOG.md updates
- Release preparation and validation
- Publishing to crates.io
- GitHub Release creation
- Backwards compatibility checks

## Versioning

Follow [Semantic Versioning 2.0.0](https://semver.org/):

- **MAJOR** (X.0.0): Breaking API changes
- **MINOR** (0.X.0): New features, backwards compatible
- **PATCH** (0.0.X): Bug fixes, backwards compatible

### Pre-1.0 Rules

While at `0.x.y`, minor versions may contain breaking changes:
- `0.1.x` → `0.2.0` may break API
- `0.1.0` → `0.1.1` should only be bug fixes

### Breaking Changes

**These require a major version bump (or minor bump pre-1.0)**:
- Removing/renaming public types, functions, methods, modules
- Changing function signatures
- Changing struct fields (public structs)
- Changing trait bounds that break existing impls
- Removing trait implementations
- Increasing MSRV
- Changing default features that affect compilation

**These are NOT breaking** (minor/patch acceptable):
- Adding new public items
- Adding new trait implementations
- Adding new optional features
- Deprecating items (not removing)
- Bug fixes matching documented behavior
- Performance improvements

## CHANGELOG Format

Use [Keep a Changelog](https://keepachangelog.com/) format:

```markdown
# Changelog

## [Unreleased]

### Added
- New feature description (#PR)

### Changed
- **BREAKING**: Changed X to Y (#PR)

### Deprecated
- `old_method()` - use `new_method()` instead (#PR)

### Removed
- Removed deprecated `ancient_method()` (#PR)

### Fixed
- Fixed bug in X (#PR)

### Security
- Updated dependency to fix CVE (#PR)

## [0.2.0] - 2026-02-15

### Added
...
```

### CHANGELOG Rules

1. **Every PR updates `[Unreleased]`** section
2. **Use present tense**: "Add feature" not "Added feature"
3. **Reference PR numbers**: `(#42)` at end of line
4. **Mark breaking changes**: Prefix with `**BREAKING**:`
5. **Be user-focused**: Describe impact, not implementation

## Release Checklist

Before releasing, verify:

```bash
# 1. All CI passes
gh pr checks --watch

# 2. Version updated in Cargo.toml
grep "^version" Cargo.toml

# 3. CHANGELOG has entry for new version
grep "## \[$(cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')\]" CHANGELOG.md

# 4. Documentation builds
cargo doc --no-deps

# 5. All tests pass
cargo test

# 6. No new warnings
cargo clippy

# 7. License check passes
cargo deny check
```

## Release Process

### Step 1: Create Release Branch

```bash
git fetch origin main
git checkout -b release/v0.2.0 origin/main
```

### Step 2: Update Version

Edit `Cargo.toml`:
```toml
version = "0.2.0"
```

Then update lockfile:
```bash
cargo check
```

### Step 3: Finalize CHANGELOG

Move `[Unreleased]` content to new version section:

```markdown
## [Unreleased]

## [0.2.0] - 2026-02-15

### Added
- (content from Unreleased)
```

### Step 4: Create Release PR

```bash
git add Cargo.toml Cargo.lock CHANGELOG.md
git commit -m "chore: release v0.2.0"
git push origin release/v0.2.0
gh pr create --title "Release v0.2.0" --body "Release preparation for v0.2.0

## Checklist
- [ ] Version in Cargo.toml matches
- [ ] CHANGELOG updated with release date
- [ ] All CI checks pass
- [ ] Documentation builds without warnings"
```

### Step 5: Merge PR

After PR is approved and merged, the release is **fully automated**:

1. **Auto-tag workflow** (`.github/workflows/release-tag.yml`) triggers on merge
   - Extracts version from branch name (`release/v0.2.0` → `v0.2.0`)
   - Verifies Cargo.toml version matches
   - Verifies CHANGELOG has release entry
   - Creates and pushes the git tag

2. **Publish workflow** (`.github/workflows/publish.yml`) triggers on tag push
   - Verifies version matches tag
   - Runs full CI suite
   - Publishes to crates.io

Monitor with:
```bash
# Watch the auto-tag workflow
gh run list --workflow=release-tag.yml
gh run watch

# Then watch the publish workflow
gh run list --workflow=publish.yml
gh run watch
```

### Step 6: Create GitHub Release (Optional)

The publish workflow handles crates.io. For a GitHub Release:

```bash
# Extract release notes from CHANGELOG
gh release create v0.2.0 \
  --title "v0.2.0" \
  --notes "See CHANGELOG.md for details"
```

Or use the CHANGELOG section as release notes.

**Note**: You can also create the GitHub Release manually from the Actions summary page after the publish workflow completes.

## Pre-Release Versions

For testing before official release:

```toml
version = "0.2.0-alpha.1"  # Early testing
version = "0.2.0-beta.1"   # Feature complete
version = "0.2.0-rc.1"     # Release candidate
```

Tag and publish:
```bash
git tag -a v0.2.0-alpha.1 -m "v0.2.0-alpha.1"
git push origin v0.2.0-alpha.1
```

Pre-releases publish to crates.io but aren't selected by default.

## Deprecation Policy

1. **Deprecate for at least one minor version before removing**
2. **Provide migration path in deprecation message**
3. **Document in CHANGELOG under `### Deprecated`**

Example:
```rust
#[deprecated(
    since = "0.2.0",
    note = "Use `Rotor::from_euler_angles()` instead"
)]
pub fn from_euler(roll: T, pitch: T, yaw: T) -> Self {
    Self::from_euler_angles(roll, pitch, yaw)
}
```

## Hotfix Process

For urgent fixes to a released version:

```bash
# Branch from the release tag
git checkout -b hotfix/v0.2.1 v0.2.0

# Make fix, update version to 0.2.1
# Update CHANGELOG with fix

# Create PR targeting main
gh pr create --title "Hotfix v0.2.1" --body "..."

# After merge, tag and release
git tag -a v0.2.1 -m "Hotfix v0.2.1"
git push origin v0.2.1
```

## Version Query Commands

```bash
# Current version in Cargo.toml
cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version'

# Latest published version
cargo search clifford --limit 1

# Check if version is published
VERSION="0.2.0"
cargo search clifford | grep -q "\"$VERSION\"" && echo "Published" || echo "Not published"
```

## Troubleshooting

### Publish Failed

1. Check CARGO_REGISTRY_TOKEN is set in repository secrets
2. Verify version isn't already published
3. Check for yanked dependencies
4. Review cargo publish --dry-run output

### Version Mismatch

If tag doesn't match Cargo.toml:
```bash
# Delete the incorrect tag
git tag -d v0.2.0
git push origin :refs/tags/v0.2.0

# Fix version, commit, and re-tag
```

### CHANGELOG Missing Entry

The publish workflow verifies CHANGELOG has an entry for the version being released. If missing:

1. Add the version section to CHANGELOG.md
2. Amend the release commit or create new commit
3. Update the tag to point to the new commit
