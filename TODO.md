# AtomPhys: Specific Development Tasks

## Critical Bugs

1. [ ] Fix `ac_stark.py` line 59-60: Define missing `tr_omega` variable in first transitions loop
   - Add `tr_omega = transition.angular_frequency` before line 42
   - Ensure proper unit handling consistent with second transitions loop

2. [ ] Remove debug print statement in `ac_stark.py` line 100
   - Delete `print(delta_E)` statement

3. [ ] Fix exception handling in `ac_stark.py` lines 63-64 and 97-98
   - Replace empty `except Exception: pass` blocks with proper logging or specific exception handling
   - Consider adding specific exception types (e.g., `ValueError`, `TypeError`)

## Code Improvements

4. [ ] Implement nuclear g-factor in `state.py` line 249
   - Replace placeholder `gI = 0` with actual nuclear g-factor calculation
   - Add lookup table for common isotopes' gI values

5. [ ] Fix transition type determination in `transition.py` line 196
   - Implement logic to determine E1/M1/E2 transition types from quantum numbers
   - Add validation against known transitions

6. [ ] Add proper mass number handling in `data_utils/nist.py` line 54
   - Implement storage and usage of mass number for isotope-specific calculations
   - Update atom initialization to use mass number when provided

7. [ ] Improve wavelength-to-RGB conversion in `plot_utils/plot.py` line 208
   - Replace current implementation with a proper matplotlib colormap
   - Add configuration options for color scaling

8. [ ] Refactor unit conversion in `ac_stark.py` lines 52-54 and 86-88
   - Remove redundant conversion to MHz magnitude and back to MHz units
   - Keep quantities as Pint objects throughout calculation

## Testing

9. [ ] Add tests for `ac_stark.py` functions
   - Create test cases with known AC Stark shift values from literature
   - Add tests for polarizability calculation
   - Test with various electric field configurations

10. [ ] Add tests for hyperfine structure calculations
    - Create test cases for specific atoms with known hyperfine constants
    - Verify calculated energy levels against reference data

11. [ ] Add tests for matrix element calculations
    - Test angular factors against reference values
    - Test radial matrix elements against literature values

12. [ ] Add integration tests for full calculation workflows
    - Test complete simulation of atom-light interaction
    - Test calculation chains (e.g., polarizability → AC Stark → energy levels)

## Documentation

13. [ ] Add detailed docstrings to all functions in `calc/ac_stark.py`
    - Include descriptions of parameters, return values, units
    - Add theoretical background and equations implemented
    - Add example usage

14. [ ] Complete documentation in `docs/user_guide/calculations/ac_stark_shifts.md`
    - Add missing theoretical background
    - Include example with realistic values
    - Add explanatory plots

15. [ ] Create step-by-step tutorial for calculating magic wavelengths
    - Use `examples/Magic wavelengths for Ca+ (J. Jiang et al 2017).ipynb` as basis
    - Add detailed explanations for each step
    - Compare results with published values

16. [ ] Add type hints to all function parameters and return values
    - Prioritize core classes: `Atom`, `State`, `Transition`
    - Add proper typing for Pint quantity objects
    - Include generics where appropriate

## Features

17. [ ] Implement caching for NIST database queries in `data_utils/nist.py`
    - Add persistent cache using `data_utils/cache.py`
    - Add cache invalidation logic
    - Add option to force refresh

18. [ ] Add visualization for Stark maps
    - Create function to calculate energy levels vs. electric field strength
    - Add plotting function with customizable appearance
    - Include examples in documentation

19. [ ] Complete Zeeman effect implementation in `ac_stark.py`
    - Uncomment and fix Zeeman shift code (lines 22, 25, 33)
    - Add proper magnetic field parameter handling
    - Test against known Zeeman splitting values

20. [ ] Add automatic calculation of polarizability spectra
    - Create function to scan over wavelengths
    - Add plotting utilities for polarizability spectra
    - Include calculation of magic wavelengths

## Package Organization

21. [ ] Consolidate notebooks from `_notebooks/` and `examples/`
    - Move all notebooks to a single `examples/` directory
    - Organize by application category
    - Update references in documentation

22. [ ] Add example configuration file for common calculation parameters
    - Create template with standard laser parameters
    - Add common atoms configuration
    - Include example unit preferences

23. [ ] Update package version to reflect recent changes
    - Update version in `pyproject.toml`
    - Add proper semantic versioning rules
    - Create CHANGELOG entry

24. [ ] Fix all ruff linting warnings
    - Run `poetry run ruff check .`
    - Address all style issues
    - Configure pre-commit hooks

## High Priority Tasks (Next 2 Weeks)

- Fix critical bugs in AC Stark calculation (tasks 1-3)
- Add tests for AC Stark and polarizability (task 9)
- Improve documentation for AC Stark calculations (tasks 13-14)
- Fix transition type determination (task 5)
- Fix ruff linting warnings (task 24)