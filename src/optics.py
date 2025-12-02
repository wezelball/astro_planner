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
        return self.focal_length_mm * (1.0 / self.focal_reducer)

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
        if object_size_deg is None:
            return 0.0

        fov_w, fov_h, _, _ = self.fov_degrees()
        min_fov = min(fov_w, fov_h)

        return min(1.0, object_size_deg / min_fov)

