class Instrument:

    def __init__(self, specs):
    
        self.name = specs['name']
        self.nPixelsX = specs['nPixelsX']                        # nr of pixels along x axis
        self.nPixelsY = specs['nPixelsY']                        # nr of pixels along y axis
        self.nOversampleX = specs['nOversampleX']
        self.nOversampleY = specs['nOversampleY']
        self.PhiX = specs['PhiX']                                # instrument FOV in x (half opening angle) in degrees
        self.PhiY = specs['PhiY']                                # instrument FOV in y (half opening angle) in degrees
        self.iFOV = specs['iFOV']                                  # pixel FOV in rad
        self.PixelSize = specs['PixelSize']                      # area of one pixel
        self.PixelFOV = specs['PixelFOV']

        self.computedQuantity = specs['computedQuantity']
