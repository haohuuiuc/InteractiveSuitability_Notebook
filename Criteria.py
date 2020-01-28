import arcpy
import math
import statistics
import matplotlib.pyplot as plt


class Criteria:
    def __init__(self, raster_obj):
        self.raster = raster_obj
        # Store values in a list
        x = []
        for r, c in raster_obj:
            x.append(raster_obj[r, c])
        # exclude nodata value, may take a while depending on raster size
        self.values = list(filter(lambda v: not math.isnan(v), x))
        self.transformed_values = self.values.copy()
        self.min_value = min(self.values)
        self.max_value = max(self.values)
        self.mean_value = statistics.mean(self.values)
        self.std_value = statistics.stdev(self.values)
        self.n_values = len(self.values)
        interv = (self.max_value - self.min_value) / (100 - 1)
        self.sample_values = [min(self.values) + i * interv for i in range(100)]
        self.transformed_sample_values = self.sample_values.copy()

    def stats(self):
        v_mean = statistics.mean(self.values)
        print('Mean: {}'.format(v_mean))
        print('Min: {}'.format(self.min_value))
        print('Max: {}'.format(self.max_value))

    def transform_stats(self):
        v_mean = statistics.mean(self.transformed_values)
        print('Mean: {}'.format(v_mean))
        print('Min: {}'.format(min(self.transformed_values)))
        print('Max: {}'.format(max(self.transformed_values)))

    def show_hist(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax1.hist(self.values, bins=n_bins)
        ax1.set_xlabel(self.raster.name)
        ax1.set_ylabel('Count')
        plt.title('Histogram of {}'.format(self.raster.name))
        plt.show()

    def transform(self, type, params):
        if type == 'unique':
            for idx, val in enumerate(self.transformed_values):
                if val in params['remap']:
                    self.transformed_values[idx] = params['remap'][val]

        if type == 'range':
            for idx, val in enumerate(self.transformed_values):
                for s,e in params['remap']:
                    if s < val <= e:
                        self.transformed_values[idx] = params['remap'][(s,e)]

        if type == 'continous':
            # RBF Small method
            if params['name'] == 'small':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = 5
                self.transformed_sample_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                                  for i in self.sample_values]
                self.transformed_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                           for i in self.values]

            # RBF Large method
            if params['name'] == 'large':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = -5
                self.transformed_sample_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                                  for i in self.sample_values]
                self.transformed_values = [1 / (1 + math.pow(i / params['mid_point'], params['spread']))
                                           for i in self.values]

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
                self.transformed_values = [n_std / (i - n_mean + n_std) if i > n_mean else 1 for i in self.values]

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
                self.transformed_values = [1 - n_std / (i - n_mean + n_std) if i > n_mean else 0 for i in self.values]

            # RBF Gaussion method
            if params['name'] == 'gaussian':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = math.log(10) * 4 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [math.exp(-params['spread'] * (i - params['mid_point']) ** 2) for i in
                                                  self.sample_values]
                self.transformed_values = [math.exp(-params['spread'] * (i - params['mid_point']) ** 2) for i in
                                           self.values]

            # RBF Near method
            if params['name'] == 'near':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = 36 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [1 / (1 + params['spread'] * math.pow(i - params['mid_point'], 2))
                                                  for i in self.sample_values]
                self.transformed_values = [1 / (1 + params['spread'] * math.pow(i - params['mid_point'], 2))
                                           for i in self.values]

            # RBF Linear method
            if params['name'] == 'linear':
                if 'min_x' not in params:
                    params['min_x'] = self.min_value
                if 'max_x' not in params:
                    params['max_x'] = self.max_value
                diff = params['max_x'] - params['min_x']
                y_sample = []
                y = []
                if diff > 0:  # positive slope
                    for v in self.sample_values:
                        if v < params['min_x']:
                            y_sample.append(0)
                        else:
                            if v > params['max_x']:
                                y_sample.append(1)
                            else:
                                y_sample.append((v - params['min_x']) / diff)
                else:
                    for v in self.values:
                        if v > params['min_x']:
                            y.append(0)
                        else:
                            if v < params['max_x']:
                                y.append(1)
                            else:
                                y.append((v - params['min_x']) / diff)
                self.transformed_sample_values = y_sample
                self.transformed_values = y

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
                if diff > 0:  # positive slope
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
                else:
                    for v in self.sample_values:
                        if v < params['max_x']:
                            y.append(1)
                        else:
                            if v < mid_p:
                                y.append((v - mid_p) / h_diff)
                            else:
                                if v > params['min_x']:
                                    y.append(1)
                                else:
                                    y.append((mid_p - v) / h_diff)
                self.transformed_sample_values = y_sample
                self.transformed_values = y

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
                self.transformed_values = [math.exp((i - params['in_shift']) * params['base_factor'])
                                           for i in self.transformed_values]

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
                self.transformed_values = [math.log((i - params['in_shift']) * params['base_factor'])
                                           for i in self.transformed_values]

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
                self.transformed_values = [math.pow(i - params['in_shift'], params['exponent']) for i in self.values]

            # RBF Logistic Growth method
            if params['name'] == 'logisticgrowth':
                if 'y_intercept_percent' not in params:
                    params['y_intercept_percent'] = 1
                c = 100
                a = c / params['y_intercept_percent'] - 1
                b = - math.log(a) / (0.5 * (self.max_value + self.min_value) - self.min_value)
                self.transformed_sample_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in
                                                  self.sample_values]
                self.transformed_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in self.values]

            # RBF Logistic Decay method
            if params['name'] == 'logisticdecay':
                if 'y_intercept_percent' not in params:
                    params['y_intercept_percent'] = 99
                c = 100
                a = c / params['y_intercept_percent'] - 1
                b = - math.log(a) / (0.5 * (self.max_value + self.min_value) - self.min_value)
                self.transformed_sample_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in
                                                  self.sample_values]
                self.transformed_values = [c / (1 + a * math.exp((i - self.min_value) * b)) for i in self.values]

            min_transformed_sample_value = min(self.transformed_sample_values)
            max_transformed_sample_value = max(self.transformed_sample_values)
            self.transformed_sample_values = [
                (v - min_transformed_sample_value) / (max_transformed_sample_value - min_transformed_sample_value) *
                (params['to_scale'] - params['from_scale']) + params['from_scale'] for v in
                self.transformed_sample_values]

        min_transformed_value = min(self.transformed_values)
        max_transformed_value = max(self.transformed_values)
        self.transformed_values = [(v - min_transformed_value) / (max_transformed_value - min_transformed_value) *
                                   (params['to_scale'] - params['from_scale']) + params['from_scale'] for v in
                                   self.transformed_values]

    def show_transform_hist(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax1.hist(self.transformed_values, bins=n_bins)
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
    c1 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\dem_24'))
    c2 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\landuse2002'))
    # c1.stats()
    # c1.show_hist()


    # RBF params
    transform_params_continous = {
        'name': 'mssmall',
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
    c1.transform('range', transform_params_range)
    c1.show_transform_hist()
    c1.transform_stats()

if __name__ == "__main__":
    main()
