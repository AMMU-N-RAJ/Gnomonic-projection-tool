"""
Orion Star Plate Coordinate Converter
Converts RA/DEC celestial coordinates of Orion constellation stars
to physical X,Y coordinates on a hexagonal star plate using gnomonic projection.
"""

import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─── Constants ────────────────────────────────────────────────────────────────
HEX_SIDE = 37.5            # mm, hexagon side length
PLATE_DIAMETER = 75.0       # mm
MARGIN_FACTOR = 0.90   # 90% of max scale (10% buffer from hex edge)
D_MAX_APERTURE = 1.0        # mm, max aperture for brightest star

INPUT_CSV = r"d:\Orion\task\stardata.csv"
OUTPUT_CSV = r"d:\Orion\task\orion_plate_output.csv"
OUTPUT_PNG = r"d:\Orion\task\orion_plate_diagram.png"

SQRT3 = math.sqrt(3)


# ─── Parsing ──────────────────────────────────────────────────────────────────

def parse_ra(ra_string):
    """Convert RA string 'HH MM SS.ss' to decimal degrees."""
    parts = ra_string.strip().split()
    hours = float(parts[0])
    minutes = float(parts[1])
    seconds = float(parts[2])
    return hours * 15.0 + minutes * (15.0 / 60.0) + seconds * (15.0 / 3600.0)


def parse_dec(dec_string):
    """Convert DEC string '+/-DD MM SS.s' to decimal degrees."""
    s = dec_string.strip()
    sign = -1.0 if s.startswith('-') else 1.0
    s = s.lstrip('+-')
    parts = s.split()
    degrees = float(parts[0])
    arcminutes = float(parts[1])
    arcseconds = float(parts[2])
    return sign * (degrees + arcminutes / 60.0 + arcseconds / 3600.0)


def read_star_data(filepath):
    """Read the star CSV file, handling BOM, trailing commas, and empty rows."""
    stars = []
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            row = [field.strip() for field in row]
            if not row or not row[0]:
                continue
            stars.append({
                'name': row[0],
                'hip_id': row[1],
                'ra_hms': row[2],
                'dec_hms': row[3],
                'mag': float(row[4]),
            })
    return stars


# ─── Projection ───────────────────────────────────────────────────────────────

def gnomonic_project(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """
    Gnomonic (tangent plane) projection of (ra, dec) onto plane tangent at (ra0, dec0).
    All inputs in degrees. Returns (x, y) in radians on the tangent plane.
    """
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    ra0 = math.radians(ra0_deg)
    dec0 = math.radians(dec0_deg)

    delta_ra = ra - ra0

    cos_c = (math.sin(dec0) * math.sin(dec) +
             math.cos(dec0) * math.cos(dec) * math.cos(delta_ra))

    x = math.cos(dec) * math.sin(delta_ra) / cos_c
    y = (math.cos(dec0) * math.sin(dec) -
         math.sin(dec0) * math.cos(dec) * math.cos(delta_ra)) / cos_c

    return x, y


# ─── Scaling & Boundary ──────────────────────────────────────────────────────

def is_inside_hexagon(x_mm, y_mm, side):
    """
    Check if (x_mm, y_mm) is inside a flat-top regular hexagon
    centered at origin with given side length.
    """
    ax = abs(x_mm)
    ay = abs(y_mm)
    if ax > side:
        return False
    if ay > side * SQRT3 / 2.0:
        return False
    if SQRT3 * ax + ay > SQRT3 * side:
        return False
    return True


def compute_scale_factor(projected_stars, hex_side, margin=0.90):
    """
    Find scale factor (mm/rad) so all stars fit within the hexagon, with margin.
    Uses binary search over scale values.
    """
    def all_inside(scale):
        for x_rad, y_rad in projected_stars:
            x_mm = x_rad * scale
            y_mm = y_rad * scale
            if not is_inside_hexagon(x_mm, y_mm, hex_side):
                return False
        return True

    lo, hi = 10.0, 1000.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if all_inside(mid):
            lo = mid
        else:
            hi = mid

    return lo * margin


# ─── Aperture ─────────────────────────────────────────────────────────────────

def compute_aperture(mag, mag_min, d_max):
    """
    Compute aperture diameter based on stellar magnitude.
    Brighter stars (lower mag) get larger holes.
    d = d_max * 10^(-0.2 * (mag - mag_min))
    """
    return d_max * (10.0 ** (-0.2 * (mag - mag_min)))


# ─── Output ───────────────────────────────────────────────────────────────────

def write_output_csv(stars, output_path):
    """Write results to CSV file."""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Star Name', 'RA (deg)', 'DEC (deg)',
                         'X (mm)', 'Y (mm)', 'Aperture (mm)', 'Inside Hexagon'])
        for s in stars:
            writer.writerow([
                s['name'],
                f"{s['ra_deg']:.4f}",
                f"{s['dec_deg']:.4f}",
                f"{s['x_mm']:.3f}",
                f"{s['y_mm']:.3f}",
                f"{s['aperture_mm']:.4f}",
                'Yes' if s['inside_hex'] else 'NO - WARNING'
            ])


def create_visualization(stars, hex_side, output_path):
    """Generate matplotlib diagram of star positions on the hexagonal plate."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_aspect('equal')

    # Plate circle (radius = hex_side)
    plate_circle = plt.Circle((0, 0), hex_side, fill=False,
                               linestyle='--', color='gray', linewidth=1)
    ax.add_patch(plate_circle)

    # Hexagon (flat-top orientation)
    s = hex_side
    hex_vertices = [
        (s, 0),
        (s / 2, s * SQRT3 / 2),
        (-s / 2, s * SQRT3 / 2),
        (-s, 0),
        (-s / 2, -s * SQRT3 / 2),
        (s / 2, -s * SQRT3 / 2),
    ]
    hex_patch = plt.Polygon(hex_vertices, closed=True, fill=False,
                             edgecolor='black', linewidth=2)
    ax.add_patch(hex_patch)

    # Plot stars
    for star in stars:
        x, y = star['x_mm'], star['y_mm']
        aperture = star['aperture_mm']

        # Scale aperture for visible marker size (points^2 for scatter)
        visual_size = (aperture * 8) ** 2
        ax.scatter(x, y, s=visual_size, c='gold', edgecolors='black',
                   linewidths=0.5, zorder=5)

        ax.annotate(star['name'], (x, y), textcoords="offset points",
                    xytext=(5, 5), fontsize=6, ha='left')

    ax.set_xlim(-80, 80)
    ax.set_ylim(-80, 80)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_title('Orion Star Plate - Hexagonal Active Area')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Diagram saved to {output_path}")
    plt.close()


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    # Step 1: Read CSV
    stars = read_star_data(INPUT_CSV)
    print(f"Loaded {len(stars)} stars from {INPUT_CSV}\n")

    # Step 2: Convert RA/DEC to decimal degrees
    for s in stars:
        s['ra_deg'] = parse_ra(s['ra_hms'])
        s['dec_deg'] = parse_dec(s['dec_hms'])

    # Step 3: Compute reference center
    ra0 = sum(s['ra_deg'] for s in stars) / len(stars)
    dec0 = sum(s['dec_deg'] for s in stars) / len(stars)
    print(f"Reference center: RA0 = {ra0:.4f} deg, DEC0 = {dec0:.4f} deg")

    # Step 4: Gnomonic projection
    for s in stars:
        s['x_proj'], s['y_proj'] = gnomonic_project(
            s['ra_deg'], s['dec_deg'], ra0, dec0)

    # Step 5: Scale to physical plate
    projected = [(s['x_proj'], s['y_proj']) for s in stars]
    scale = compute_scale_factor(projected, HEX_SIDE, MARGIN_FACTOR)
    print(f"Scale factor: {scale:.2f} mm/rad\n")

    for s in stars:
        s['x_mm'] = s['x_proj'] * scale
        s['y_mm'] = s['y_proj'] * scale

    # Step 6: Boundary check
    for s in stars:
        s['inside_hex'] = is_inside_hexagon(s['x_mm'], s['y_mm'], HEX_SIDE)
        if not s['inside_hex']:
            print(f"WARNING: {s['name']} at ({s['x_mm']:.3f}, {s['y_mm']:.3f}) "
                  f"is OUTSIDE the hexagonal boundary!")

    # Step 7: Compute apertures
    mag_min = min(s['mag'] for s in stars)
    for s in stars:
        s['aperture_mm'] = compute_aperture(s['mag'], mag_min, D_MAX_APERTURE)

    # Step 8: Print results table
    print(f"{'Star':<30} {'RA(deg)':>9} {'DEC(deg)':>9} "
          f"{'X(mm)':>8} {'Y(mm)':>8} {'Aper(mm)':>9} {'InHex':>6}")
    print("-" * 90)
    for s in stars:
        print(f"{s['name']:<30} {s['ra_deg']:9.4f} {s['dec_deg']:9.4f} "
              f"{s['x_mm']:8.3f} {s['y_mm']:8.3f} {s['aperture_mm']:9.4f} "
              f"{'Yes' if s['inside_hex'] else 'NO':>6}")

    # Step 9: Write output CSV
    write_output_csv(stars, OUTPUT_CSV)
    print(f"\nResults written to {OUTPUT_CSV}")

    # Step 10: Generate visualization
    create_visualization(stars, HEX_SIDE, OUTPUT_PNG)


if __name__ == '__main__':
    main()
