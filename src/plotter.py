import numpy as np

def build_ir_spectrum(
    conformers,
    freq_min=400.0,
    freq_max=4000.0,
    resolution=1.0,
    fwhm=15.0,
    scale_factor=0.97,
):
    """
    Build a Boltzmann-weighted IR spectrum using Gaussian broadening.

    Args:
        conformers: list of conformer dicts with keys:
            - ir_freqs (cm^-1)
            - ir_intensities (km/mol)
            - weight
        freq_min, freq_max: spectral range (cm^-1)
        resolution: grid spacing (cm^-1)
        fwhm: Gaussian full width at half max (cm^-1)
        scale_factor: harmonic frequency scaling factor

    Returns:
        freqs: np.ndarray
        spectrum: np.ndarray
    """
    if not all(np.isfinite(value) for value in (freq_min, freq_max, resolution, fwhm, scale_factor)):
        raise ValueError("Spectrum settings must be finite")
    if freq_min >= freq_max or resolution <= 0 or fwhm <= 0 or scale_factor <= 0:
        raise ValueError("Invalid spectrum range, resolution, FWHM, or scale factor")
    freqs = np.arange(freq_min, freq_max + resolution * 0.5, resolution)
    spectrum = np.zeros_like(freqs, dtype=float)

    # Convert FWHM → sigma
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    for c in conformers:
        if "ir_freqs" not in c or "ir_intensities" not in c:
            continue

        weight = c.get("weight", 0.0)
        if weight <= 0.0:
            continue

        for nu, intensity in zip(c["ir_freqs"], c["ir_intensities"]):
            nu_scaled = nu * scale_factor
            if not np.isfinite(nu) or not np.isfinite(intensity) or intensity < 0:
                continue
            gaussian = intensity * np.exp(
                -0.5 * ((freqs - nu_scaled) / sigma) ** 2
            ) / (sigma * np.sqrt(2.0 * np.pi))
            spectrum += weight * gaussian

    return freqs, spectrum
