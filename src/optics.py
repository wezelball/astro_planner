# src/optics.py
import math

class Optics:
    """
    Represents telescope + camera system.
    Computes FOV, pixel scale, and FOV fill fraction.
    """

    def __init__(
        self,
        telescope_aperture_mm,
        telescope_focal_length_mm,
        camera_sensor_width_mm,
        camera_sensor_height_mm,
        camera_pixel_width,
        camera_pixel_height,
        focal_reducer=1.0
    ):
        self.aperture_mm = telescope_aperture_mm
        self.focal_length_mm = telescope_focal_length_mm
        self.sensor_width_mm = camera_sensor_width_mm
        self.sensor_height_mm = camera_sensor_height_mm
        self.pixel_width = camera_pixel_width
        self.pixel_height = camera_pixel_height
        self.focal_reducer = focal_reducer

    # ---------------------------------------------------------
    # Effective focal length
    # ---------------------------------------------------------
    @property
    def effective_focal_length(self):
        # assume self.focal_length_mm already stores the effective focal length (applied reducer)
        return self.focal_length_mm
    
    # ---------------------------------------------------------
    # TODO - ADD DESCRIPTION
    # ---------------------------------------------------------
    def fov_diagonal_deg(self):
        """Return diagonal FOV in degrees."""
        fov_w, fov_h, _, _ = self.fov_degrees()
        return math.sqrt(fov_w**2 + fov_h**2)

    # ---------------------------------------------------------
    # TODO - ADD DESCRIPTION
    # ---------------------------------------------------------
    def pixel_scale_arcsec(self):
        """
        Return pixel scale in arcseconds/pixel: (x_arcsec, y_arcsec)
        Uses sensor size in mm -> convert pixel size to microns implicitly via *1000 factor.
        """
        # pixel_size_um = (sensor_mm * 1000.0) / pixels
        scale_x = 206.265 * (self.sensor_width_mm * 1000.0 / self.pixel_width) / self.effective_focal_length_mm
        scale_y = 206.265 * (self.sensor_height_mm * 1000.0 / self.pixel_height) / self.effective_focal_length_mm
        return (scale_x, scale_y)


    # ---------------------------------------------------------
    # FOV in degrees and arcseconds per pixel
    # ---------------------------------------------------------
    def fov_degrees(self):
        ef = self.effective_focal_length

        # FOV in degrees
        fov_w = (57.2957795 * self.sensor_width_mm / ef)
        fov_h = (57.2957795 * self.sensor_height_mm / ef)

        # Pixel scale in arcsec/pixel
        scale_w = (206.265 * (self.sensor_width_mm / self.pixel_width)) / ef
        scale_h = (206.265 * (self.sensor_height_mm / self.pixel_height)) / ef

        return fov_w, fov_h, scale_w, scale_h

    # ---------------------------------------------------------
    # Fraction of FOV the object occupies (based on major axis)
    # ---------------------------------------------------------
    def fov_fill_fraction(self, object_size_deg):
        """
        Strict definition: fraction = object_major_axis_deg / diagonal_fov_deg.
        object_size_deg must be a float (degrees). If object_size_deg is None, caller should decide (here we raise).
        """
        print(f"ojject_size_deg: {object_size_deg}")
        if object_size_deg is None:
            # We raise here to make caller intentionally handle missing sizes.
            raise ValueError("object_size_deg is None; Optics.fov_fill_fraction requires a numeric size in degrees")

        diag = self.fov_diagonal_deg()

        if diag <= 0:
            return 0.0
        return object_size_deg / diag

