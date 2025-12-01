import math

# ----------------------------
# TELESCOPE
# ----------------------------
class Telescope:
    def __init__(self, name, aperture_mm, focal_length_mm):
        self.name = name
        self.aperture_mm = aperture_mm
        self.focal_length_mm = focal_length_mm


class Camera:
    def __init__(self, name, sensor_width_mm, sensor_height_mm, pixel_width, pixel_height):
        self.name = name
        self.sensor_width_mm = sensor_width_mm
        self.sensor_height_mm = sensor_height_mm
        self.pixel_width = pixel_width
        self.pixel_height = pixel_height


class Optics:
    def __init__(self, telescope: Telescope, camera: Camera, focal_reducer=1.0):
        self.telescope = telescope
        self.camera = camera
        self.focal_reducer = focal_reducer  # single optional reducer

    @property
    def effective_focal_length_mm(self):
        return self.telescope.focal_length_mm * self.focal_reducer

    @property
    def fov_deg(self):
        fov_x = 2 * math.degrees(math.atan(self.camera.sensor_width_mm / (2 * self.effective_focal_length_mm)))
        fov_y = 2 * math.degrees(math.atan(self.camera.sensor_height_mm / (2 * self.effective_focal_length_mm)))
        return fov_x, fov_y

    @property
    def pixel_scale_arcsec(self):
        scale_x = 206.265 * (self.camera.sensor_width_mm * 1000 / self.camera.pixel_width) / self.effective_focal_length_mm
        scale_y = 206.265 * (self.camera.sensor_height_mm * 1000 / self.camera.pixel_height) / self.effective_focal_length_mm
        return scale_x, scale_y

    def fov_fill_fraction(self, object_size_deg):
        fov_x, fov_y = self.fov_deg
        fov_max = max(fov_x, fov_y)
        return object_size_deg / fov_max

