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

    def show_hist(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax1.hist(self.values, bins=n_bins)
        ax1.set_xlabel(self.raster.name)
        ax1.set_ylabel('Count')
        plt.title('Histogram of {}'.format(self.raster.name))
        plt.show()

    def transform(self, type, params):
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
                self.transformed_sample_values = [n_std/ (i - n_mean + n_std) if i > n_mean else 1 for i in self.sample_values]
                self.transformed_values = [n_std/ (i - n_mean + n_std) if i > n_mean else 1 for i in self.values]

            # RBF MSLarge method
            if params['name'] == 'mslarge':
                if 'mean_multiplier' not in params:
                    params['mean_multiplier'] = 1
                if 'std_multiplier' not in params:
                    params['std_multiplier'] = 1
                n_mean = params['mean_multiplier'] * self.mean_value
                n_std = params['std_multiplier'] * self.std_value
                self.transformed_sample_values = [1-n_std/ (i - n_mean + n_std) if i > n_mean else 0 for i in self.sample_values]
                self.transformed_values = [1-n_std/ (i - n_mean + n_std) if i > n_mean else 0 for i in self.values]

            # RBF Gaussion method
            if params['name'] == 'gaussian':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = math.log(10) * 4 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [math.exp(-params['spread'] * (i - params['mid_point']) ** 2) for i in self.sample_values]
                self.transformed_values = [math.exp(-params['spread'] * (i - params['mid_point']) ** 2) for i in self.values]

            # RBF Near method
            if params['name'] == 'near':
                if 'mid_point' not in params:
                    params['mid_point'] = (self.max_value + self.min_value) / 2
                if 'spread' not in params:
                    params['spread'] = 36 / math.pow(params['mid_point'] - self.min_value, 2)
                self.transformed_sample_values = [ 1 / ( 1+ params['spread'] * math.pow(i - params['mid_point'], 2))
                                                   for i in self.sample_values]
                self.transformed_values = [ 1 / ( 1+ params['spread'] * math.pow(i - params['mid_point'], 2))
                                                   for i in self.values]




    def show_transform_plot(self, scale_min=1, scale_max=10, n_bins=20):
        self.transformed_sample_values = [(v - self.min_value) / (self.max_value - self.min_value) *
                                          (scale_max - scale_min) + scale_min for v in self.transformed_sample_values]

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
    # c1.stats()
    # c1.show_hist()

    # RBF small / large / gaussian / near
    transform_params = {
        'name': 'near'
    }

    c1.transform('continous', transform_params)
    c1.show_transform_plot()


if __name__ == "__main__":
    main()
