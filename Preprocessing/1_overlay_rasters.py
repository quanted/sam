import os
import numpy as np
import gdal
import math
import pickle
from Preprocessing.utilities import nhd_states
import time

from Tool.functions import MemoryMatrix


class Envelope(object):
    """ Object representing a simple bounding rectangle, used primarily to measure raster overlap """

    def __init__(self, left, right, bottom, top):
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

    #  Returns the rectangle corresponding to the overlap with another Envelope object
    def overlap(self, r2):
        def range_overlap(a_min, a_max, b_min, b_max):
            return (a_min <= b_max) and (b_min <= a_max)

        if not all((range_overlap(self.left, self.right, r2.left, r2.right),
                    range_overlap(self.bottom, self.top, r2.bottom, r2.top))):
            return None
        else:
            left, right = sorted([self.left, self.right, r2.left, r2.right])[1:3]
            bottom, top = sorted([self.bottom, self.top, r2.bottom, r2.top])[1:3]
        return Envelope(left, right, bottom, top)

    @property
    def area(self):
        return abs(self.top - self.bottom) * abs(self.right - self.left)

    def __repr__(self):
        return "Rectangle(left: {}, right: {}, top: {}, bottom: {}".format(self.left, self.right, self.top, self.bottom)

    def __eq__(self, other):
        if (self.left, self.right, self.bottom, self.top) == (other.left, other.right, other.bottom, other.top):
            return True
        else:
            return False


class Raster(object):
    """ Wrapper for an ESRI raster grid and function for reading into an array """

    def __init__(self, path, no_data=None, alias=None, no_data_to=0):
        self.path = path
        self.no_data = no_data
        self.obj = gdal.Open(path)
        self.gt = self.obj.GetGeoTransform()
        self.cell_size = int(self.gt[1])
        self.tl = (self.gt[0], self.gt[3])  # top left
        self.size = (self.obj.RasterXSize * self.cell_size, self.obj.RasterYSize * self.cell_size)
        self.shape = Envelope(self.tl[0], self.tl[0] + self.size[0], self.tl[1] - self.size[1], self.tl[1])
        self.obj.GetRasterBand(1).SetNoDataValue(0.0)

        self.envelope = None
        self._array = None
        self.alias = alias

    def array(self, envelope, zero_min=True, dtype=np.int64):
        if self._array is None or (envelope != self.envelope):
            offset_x = (envelope.left - self.shape.left)
            offset_y = (self.shape.top - envelope.top)
            x_max = (envelope.right - envelope.left)
            y_max = (envelope.top - envelope.bottom)
            bounds = map(lambda x: int(x / self.cell_size), (offset_x, offset_y, x_max, y_max))
            self._array = self.obj.ReadAsArray(*bounds)
            if self.alias:
                self._array = np.array([self.alias.get(val, 0) for val in self._array.flat]).reshape(self._array.shape)
            if zero_min:
                self._array[self._array < 0] = 0
            if self.no_data:
                self._array[self._array == self.no_data] = 0
            self.envelope = envelope
        return dtype(self._array)

    def build_alias(self, tile_size=500000, pickle_file=None, overwrite=False):
        if overwrite or not os.path.exists(pickle_file):
            print('Building an alias for {}...'.format(self.path))
            all_values = set()
            tiles = make_tiles(self.shape, tile_size)
            for i, tile in enumerate(tiles):
                print(i, len(tiles))
                all_values |= set(np.unique(self.array(tile)))
            self.alias = dict(zip(sorted(all_values), np.arange(len(all_values))))
            with open(pickle_file, 'wb') as f:
                pickle.dump(self.alias, f)
        else:
            with open(pickle_file, 'rb') as f:
                self.alias = pickle.load(f)

    @property
    def max_val(self):
        if self.alias is None:
            return int(self.obj.GetRasterBand(1).GetMaximum())
        else:
            return max(self.alias.values())

    @property
    def precision(self):
        return 10 ** int(math.ceil(math.log10(self.max_val)))


# Divides a bounding envelope into smaller tiles for processing on less powerful computers
def make_tiles(envelope, tile_size):
    if tile_size == 'max':
        return [envelope]
    else:
        h = list(range(int(envelope.left), int(envelope.right), tile_size)) + [envelope.right]
        v = list(range(int(envelope.bottom), int(envelope.top), tile_size)) + [envelope.top]
        return [Envelope(h[i], h[i + 1], v[j], v[j + 1]) for i in range(len(h) - 1) for j in range(len(v) - 1)]


def allocate(all_rasters, tile='max', cell_size=30):
    """ Allocates raster classes to a set of overlapping zones """

    def overlay_tile(rasters):
        combined_array = None
        for i, raster in enumerate(rasters):

            # Pull the array for the tile
            local_zone = raster.array(tile)
            if not local_zone.any():  # If there is no data for the raster in this tile, pull the plug
                return

            # If the raster has a larger cell size than the specified cell size, adjust accordingly
            cellsize_adjust = raster.cell_size / cell_size
            if cellsize_adjust > 1:
                local_zone = local_zone.repeat(cellsize_adjust, axis=0).repeat(cellsize_adjust, axis=1)

            # Add the adjusted zone raster to the combined array.

            if combined_array is None:
                combined_array = (local_zone * raster.adjust)
            else:
                try:
                    combined_array += (local_zone * raster.adjust)
                except ValueError:
                    a0_min = min((local_zone.shape[0], combined_array.shape[0]))
                    a1_min = min((local_zone.shape[1], combined_array.shape[1]))
                    combined_array[:a0_min, :a1_min] += (local_zone[:a0_min, :a1_min] * raster.adjust)
                    print("Shape mismatch: {} vs {}. Might be some missing pixels".format(combined_array.shape,
                                                                                          (a0_min, a1_min)))

        return combined_array

    # Overlap rasters and create envelope covering common areas
    overlap_area = None
    for i, raster in enumerate(all_rasters):
        overlap_area = raster.shape if overlap_area is None else raster.shape.overlap(overlap_area)
        assert overlap_area, "Zone and allocation rasters do not overlap"

    # Divide the overlap area into tiles to aid in processing
    tiles = make_tiles(overlap_area, tile)

    # Sort by precision and rearrange index
    old_index = [raster.path for raster in all_rasters]
    all_rasters = sorted(all_rasters, key=lambda x: x.precision, reverse=True)
    new_index = [old_index.index(raster.path) for raster in all_rasters]
    header = np.array(["weather_grid", "cdl", "comid", "mukey"])[new_index]

    # Multiply the zone raster by an adjustment factor based on precision
    for i, raster in enumerate(all_rasters):
        raster.adjust = np.prod([r.precision for r in all_rasters[i + 1:]], dtype=np.int64)

    # Iterate through tiles
    finished = None
    for j, tile in enumerate(tiles):

        print("\tAllocating ({}/{})".format(j + 1, len(tiles)))
        tile_array = overlay_tile(all_rasters)

        if tile_array is not None:
            # Break down combined array to get zones and classes
            values, counts = np.unique(tile_array.flat, return_counts=True)  # Cells in each zone and class
            zones = np.zeros((len(all_rasters), counts.size), dtype=np.int64)

            # Extract values from long combined integes
            for i, r in enumerate(all_rasters):
                zones[i] = np.int64(values / r.adjust)
                values = values - np.uint64(zones[i] * r.adjust)

            counts *= (cell_size ** 2)  # Convert cell count to area

            # Filter out combinations with no cells
            final = np.vstack((zones, counts))[:, (np.prod(zones, axis=0) > 0)]

            # Convert back alias
            for i, raster in enumerate(all_rasters):
                if raster.alias is not None:
                    convert_back = np.array(list(zip(*sorted(raster.alias.items(), key=lambda x: x[1])))[0])
                    final[i] = convert_back[final[i]]

            # Append to running
            if finished is not None:
                finished = np.hstack((finished, final))
            else:
                finished = final

    return finished, header


def save_output(out_table, out_file, header):
    np.savez_compressed(out_file, data=out_table, columns=header)


def main():
    # Designating a fixed dir here for tool development
    root_dir = os.path.join(r"C:\\", "Users", "shelly", "Documents", "SAM", "Trip", "Output")
    nhd_path = os.path.join(root_dir, "NHDCatchments", "region{}")
    weather_path = os.path.join(root_dir, "ProjectedLayers", "weather_grid")
    soil_path = os.path.join(root_dir, "SSURGO", "{}")
    cdl_path = os.path.join(root_dir, "ProjectedLayers", "cdl{}")
    out_file = os.path.join("..", "bin", "Preprocessed", "Combos", "{}_{}")

    # Create output directory
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    years = range(2010, 2016)
    regions = ['07'] + sorted(nhd_states.keys())
    overwrite = False

    weather_raster = Raster(weather_path)

    for year in years:
        cdl_raster = Raster(cdl_path.format(year))

        for region in regions:
            print("Processing Region {}...".format(region))
            region_alias = os.path.join('..', 'Development', 'aliases', 'nhd{}.p'.format(region))
            nhd_raster = Raster(nhd_path.format(region))
            nhd_raster.build_alias(pickle_file=region_alias)
            out_table = None
            header = None

            if overwrite or not os.path.exists(out_file.format(region, year)):
                for state in sorted(nhd_states[region]):
                    print("\tProcessing {}...".format(state))
                    start = time.time()
                    state_table = os.path.join('..', "Development", "progress",
                                               "{}_{}_{}.npy".format(state, region, year))
                    state_alias = os.path.join('..', 'Development', 'aliases', 'soils{}.p'.format(state))
                    if overwrite or not os.path.exists(state_table):
                        start = time.time()
                        soil_raster = Raster(soil_path.format(state), no_data=-2147483647)
                        soil_raster.build_alias(pickle_file=state_alias)
                        new_table, header = allocate([weather_raster, cdl_raster, nhd_raster, soil_raster], tile=300000)
                        np.save(state_table, new_table, header)
                    else:
                        new_table = np.load(state_table)

                    if out_table is None:
                        out_table = new_table.T
                    else:
                        out_table = np.concatenate([out_table, new_table.T], axis=0)
                    print("\tFinished in {} minutes".format(int((time.time() - start) / 60)))
                if header is None:
                    header = ["comid", "weather_grid", "mukey", "cdl"]
                header += ['area']
                save_output(out_table, out_file.format(region, year), header)


main()
