import arcpy
import math
import matplotlib.pyplot as plt


class Criteria:
    def __init__(self, raster_obj):
        self.raster = raster_obj
        raster_info = raster_obj.getRasterInfo()
        raster_info.setPixelType('F32')
        self.transformed_raster = arcpy.Raster(raster_info)
        self.scaled_transformed_raster = arcpy.Raster(raster_info)
        self.min_value = raster_obj.minimum
        self.max_value = raster_obj.maximum
        self.mean_value = raster_obj.mean
        self.std_value = raster_obj.standardDeviation
        self.values = self.get_raster_values(raster_obj)
        interv = (max(self.values) - min(self.values)) / (100 - 1)
        self.sample_values = [min(self.values) + i * interv for i in range(100)]
        self.transformed_sample_values = self.sample_values.copy()

    def get_raster_values(self, raster):
        x = []
        for r, c in raster:
            x.append(raster[r, c])
        # exclude nodata value, may take a while depending on raster size
        return list(filter(lambda v: not math.isnan(v), x))

    def show_stats(self, raster):
        print('Mean: {}'.format(raster.mean))
        print('Min: {}'.format(raster.minimum))
        print('Max: {}'.format(raster.maximum))

    def show_hist(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax1.hist(self.values, bins=n_bins)
        ax1.set_xlabel(self.raster.name)
        ax1.set_ylabel('Count')
        plt.title('Histogram of {}'.format(self.raster.name))
        plt.show()

    def transform(self, type, params):
        if type == 'unique':
            with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                for i, j in rci:
                    old_cell_value = self.raster[i, j]
                    if old_cell_value in params['remap']:
                        self.transformed_raster[i, j] = params['remap'][old_cell_value]

            for idx, val in enumerate(self.transformed_values):
                if val in params['remap']:
                    self.transformed_values[idx] = params['remap'][val]

        if type == 'range':
            with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                for i, j in rci:
                    old_cell_value = self.raster[i, j]
                    for s, e in params['remap']:
                        if s < old_cell_value <= e:
                            self.transformed_raster[i, j] = params['remap'][(s, e)]

        if type == 'continous':
            # RBF Small method
            if params['name'] == 'small':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = 5
                self.transformed_sample_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                                  for i in self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i,j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = 1 / (1 + math.pow(old_cell_value / params['mid_point'], params['spread']))

            # RBF Large method
            if params['name'] == 'large':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = -5
                self.transformed_sample_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                                  for i in self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = 1 / (
                                    1 + math.pow(old_cell_value / params['mid_point'], params['spread']))

            # RBF MSSmall method
            if params['name'] == 'mssmall':
                if 'mean_multiplier' not in params:
                    params['mean_multiplier'] = 1
                if 'std_multiplier' not in params:
                    params['std_multiplier'] = 1
                n_mean = params['mean_multiplier'] * self.mean_value
                n_std = params['std_multiplier'] * self.std_value
                self.transformed_sample_values = [n_std / (i - n_mean + n_std) if i > n_mean else 1 for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i,j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = n_std / (old_cell_value - n_mean + n_std) if old_cell_value > n_mean else 1

            # RBF MSLarge method
            if params['name'] == 'mslarge':
                if 'mean_multiplier' not in params:
                    params['mean_multiplier'] = 1
                if 'std_multiplier' not in params:
                    params['std_multiplier'] = 1
                n_mean = params['mean_multiplier'] * self.mean_value
                n_std = params['std_multiplier'] * self.std_value
                self.transformed_sample_values = [1 - n_std / (i - n_mean + n_std) if i > n_mean else 0 for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = old_cell_value - n_std / (old_cell_value - n_mean + n_std) if old_cell_value > n_mean else 0

            # RBF Gaussion method
            if params['name'] == 'gaussian':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = math.log(10) * 4 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [math.exp(-params['spread'] * (i - params['mid_point']) ** 2) for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = math.exp(-params['spread'] * (old_cell_value - params['mid_point']) ** 2)

            # RBF Near method
            if params['name'] == 'near':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = 36 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [1 / (1 + params['spread'] * math.pow(i - params['mid_point'], 2))
                                                  for i in self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = 1 / (1 + params['spread'] * math.pow(old_cell_value - params['mid_point'], 2))

            # RBF Linear method
            if params['name'] == 'linear':
                if 'min_x' not in params:
                    params['min_x'] = self.min_value
                if 'max_x' not in params:
                    params['max_x'] = self.max_value
                diff = params['max_x'] - params['min_x']
                y_sample = []
                y = []
                # positive slope, assume always the case
                for v in self.sample_values:
                    if v < params['min_x']:
                        y_sample.append(0)
                    else:
                        if v > params['max_x']:
                            y_sample.append(1)
                        else:
                            y_sample.append((v - params['min_x']) / diff)
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        if old_cell_value < params['min_x']:
                            self.transformed_raster[i, j] = 0
                        else:
                            if old_cell_value > params['max_x']:
                                self.transformed_raster[i, j] = 1
                            else:
                                self.transformed_raster[i, j] = (old_cell_value - params['min_x']) / diff
                self.transformed_sample_values = y_sample

            # RBF Symmetric Linear method
            if params['name'] == 'symmetriclinear':
                if 'min_x' not in params:
                    params['min_x'] = self.min_value
                if 'max_x' not in params:
                    params['max_x'] = self.max_value
                diff = params['max_x'] - params['min_x']
                h_diff = 0.5 * diff
                mid_p = params['min_x'] + h_diff
                y_sample = []
                y = []
                for v in self.sample_values:
                    if v < params['min_x']:
                        y_sample.append(0)
                    else:
                        if v < mid_p:
                            y_sample.append((v - params['min_x']) / h_diff)
                        else:
                            if v > params['max_x']:
                                y_sample.append(0)
                            else:
                                y_sample.append((params['max_x'] - v) / h_diff)
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        if old_cell_value < params['min_x']:
                            self.transformed_raster[i, j] = 0
                        else:
                            if old_cell_value < mid_p:
                                self.transformed_raster[i, j] = (old_cell_value - params['min_x']) / h_diff
                            else:
                                if old_cell_value > params['max_x']:
                                    self.transformed_raster[i, j] = 0
                                else:
                                    self.transformed_raster[i, j] = (params['max_x'] - old_cell_value) / h_diff
                self.transformed_sample_values = y_sample

            # RBF Exponential method
            if params['name'] == 'exponential':
                if 'in_shift' not in params:
                    params['in_shift'] = (self.min_value * math.log(params['to_scale']) - self.max_value *
                                          math.log(params['from_scale'])) / (math.log(params['to_scale']) -
                                                                             math.log(params['from_scale']))
                if 'base_factor' not in params:
                    params['base_factor'] = (math.log(params['to_scale']) - math.log(params['from_scale'])) / \
                                            (self.max_value - self.min_value)
                self.transformed_sample_values = [math.exp((i - params['in_shift']) * params['base_factor'])
                                                  for i in self.transformed_sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = math.exp((old_cell_value - params['in_shift']) * params['base_factor'])

            # RBF Logarithm method
            if params['name'] == 'logarithm':
                if 'in_shift' not in params:
                    params['in_shift'] = (self.min_value * math.exp(params['to_scale']) - self.max_value *
                                          math.exp(params['from_scale'])) / (math.exp(params['to_scale']) -
                                                                             math.exp(params['from_scale']))
                if 'base_factor' not in params:
                    params['base_factor'] = (math.exp(params['to_scale']) - math.exp(params['from_scale'])) / \
                                            (self.max_value - self.min_value)
                self.transformed_sample_values = [math.log((i - params['in_shift']) * params['base_factor'])
                                                  for i in self.transformed_sample_values]

                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = math.log((old_cell_value - params['in_shift']) * params['base_factor'])

            # RBF Power method
            if params['name'] == 'power':
                if params['from_scale'] == 0:
                    in_shift = self.min_value
                    if params['to_scale'] <= 1:
                        exponent = 1
                    else:
                        exponent = math.log(params['to_scale']) / (self.max_value - in_shift)
                elif params['from_scale'] == 1:
                    in_shift = self.min_value - 1
                    exponent = math.log(params['to_scale']) / math.log(self.max_value - in_shift)
                else:
                    in_shift = self.min_value
                    exponent = 2
                if 'in_shift' not in params:
                    params['in_shift'] = in_shift
                if 'exponent' not in params:
                    params['exponent'] = exponent

                self.transformed_sample_values = [math.pow(i - params['in_shift'], params['exponent']) for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = math.pow(old_cell_value - params['in_shift'], params['exponent'])

            # RBF Logistic Growth method
            if params['name'] == 'logisticgrowth':
                if 'y_intercept_percent' not in params:
                    params['y_intercept_percent'] = 1
                c = 100
                a = c / params['y_intercept_percent'] - 1
                b = - math.log(a) / (0.5 * (self.max_value + self.min_value) - self.min_value)
                self.transformed_sample_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = c / (1 + a * math.exp((old_cell_value - self.min_value) * b))

            # RBF Logistic Decay method
            if params['name'] == 'logisticdecay':
                if 'y_intercept_percent' not in params:
                    params['y_intercept_percent'] = 99
                c = 100
                a = c / params['y_intercept_percent'] - 1
                b = - math.log(a) / (0.5 * (self.max_value + self.min_value) - self.min_value)
                self.transformed_sample_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in
                                                  self.sample_values]
                with arcpy.sa.RasterCellIterator({'rasters': [self.raster, self.transformed_raster]}) as rci:
                    for i, j in rci:
                        old_cell_value = self.raster[i, j]
                        self.transformed_raster[i, j] = c / (1 + a * math.exp((old_cell_value - self.min_value) * b))


            arcpy.CalculateStatistics_management(self.transformed_raster)
            min_transformed_value = self.transformed_raster.minimum
            max_transformed_value = self.transformed_raster.maximum
            self.transformed_sample_values = [
                (v - min_transformed_value) / (max_transformed_value - min_transformed_value) *
                (params['to_scale'] - params['from_scale']) + params['from_scale'] for v in
                self.transformed_sample_values]

            # Final recale
            with arcpy.sa.RasterCellIterator({'rasters': [self.scaled_transformed_raster, self.transformed_raster]}) as rci:
                for r, c in rci:
                    self.scaled_transformed_raster[r, c] = (self.transformed_raster[r, c] - min_transformed_value) / (
                        max_transformed_value - min_transformed_value) * (params['to_scale'] - params['from_scale']) + \
                                            params['from_scale']
            self.transformed_raster = self.scaled_transformed_raster

        # Calculate statistics
        arcpy.CalculateStatistics_management(self.transformed_raster)

    def show_transformed_hist(self, n_bins=20):
        values = self.get_raster_values(self.transformed_raster)
        fig, ax1 = plt.subplots()
        ax1.hist(values, bins=n_bins)
        ax1.set_xlabel(self.raster.name)
        ax1.set_ylabel('Count')
        plt.title('Histogram of transformed {}'.format(self.raster.name))
        plt.show()

    def show_transform_plot(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.hist(self.values, bins=n_bins)
        ax2.plot(self.sample_values, self.transformed_sample_values, color='r')

        ax1.set_xlabel('X data')
        ax2.set_ylabel('Transformed value')
        ax1.set_ylabel('Count')
        plt.title('Transformed plot of {}'.format(self.raster.name))
        plt.show()


def main():
    # test code here
    c1 = Criteria(arcpy.Raster(r'\\haohu\share\\SuitabilityData\\dem_24'))
    c2 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\landuse2002'))
    # c1.stats()
    # c1.show_hist()

    # RBF params
    transform_params_continous = {
        'name': 'logisticdecay',
        'from_scale': 1,
        'to_scale': 10
    }

    # Unique transform params
    transform_params_unique = {
        'from_scale': 1,
        'to_scale': 10,
        'remap': {
            5: 1,
            7: 1,
            11: 1,
            12: 1,
            13: 1,
            14: 1,
            17: 1,
            24: 1,
            41: 10,
            42: 10,
            43: 10,
            61: 1,
            62: 1,
            211: 1,
            212: 1
        }
    }
    # Range transform params
    transform_params_range = {
        'from_scale': 1,
        'to_scale': 10,
        'remap': {
            (587, 935.8): 1,
            (935.8, 1283.6): 2,
            (1283.6, 1631.4): 3,
            (1631.4, 1979.2): 4,
            (1979.2, 2327): 5,
            (2327, 2674.8): 6,
            (2674.8, 3022.6): 7,
            (3022.6, 3370.4): 8,
            (3370.4, 3718.2): 9,
            (3718.2, 4066): 10,
        }
    }
    # c1.transform('continous', transform_params_continuous)
    # c1.show_transform_plot()
    # c1.transform_stats()
    # c2.transform('unique', transform_params_unique)
    # c2.show_transform_hist()
    # c2.transform_stats()
    # c1.transform('range', transform_params_range)
    c1.transform('continous', transform_params_continous)
    # c1.show_hist()
    # c1.show_transformed_hist()
    c1.show_transform_plot()
    c1.show_stats(c1.transformed_raster)


if __name__ == "__main__":
    main()
