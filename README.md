# Orion Star Plate Coordinate Converter

This project converts Orion constellation star RA/DEC coordinates into X/Y positions on a hexagonal star plate using a gnomonic projection. It also computes drill aperture sizes based on stellar magnitude and can generate a plot of the plate.

## Files

- orion.py: Main script for converting coordinates, generating the output CSV, and rendering a PNG diagram.
- orion2.py: Alternate version with different plate dimensions.
- stardata.csv: Input star catalog (RA/DEC in HMS/DMS, magnitude).
- orion_plate_output.csv: Example output produced by the script.

## Requirements

- Python 3.8+
- matplotlib

Install dependencies:

```
pip install matplotlib
```

## Usage

Run the main script from the project folder:

```
python orion.py
```

Outputs:

- orion_plate_output.csv: Computed positions and aperture sizes.
- orion_plate_diagram.png: Visualization of the hexagonal plate and star locations.

## Notes

- Input and output paths are currently absolute Windows paths in the scripts. If you move the project, update the `INPUT_CSV`, `OUTPUT_CSV`, and `OUTPUT_PNG` constants.
- The plate is modeled as a flat-top hexagon and all units are millimeters.
