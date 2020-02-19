from Criteria import *


class SuitabilityModel:
    def __init__(self, weight_method='multiplier', from_scale=1, to_scale=10):
        self.criteria = []
        self.weight = []
        self.weight_method = weight_method
        self.suitability_values = []
        self.from_scale = from_scale
        self.to_scale = to_scale
        self.suitability_map = None
        arcpy.CheckOutExtension("Spatial")

    def add_criteria(self, criterion, weight=1):
        self.criteria.append(criterion)
        self.weight.append(weight)

    def calculate(self):
        raster_weight_list = []

        for i in range(len(self.criteria)):
            if self.weight_method == 'multiplier':
                raster_weight_list.append([self.criteria[i].transformed_raster, 'VALUE', self.weight[i]])
            else:
                raster_weight_list.append([self.criteria[i].transformed_raster, 'VALUE', self.weight[i]/ 100])
        WSumTableObj  = arcpy.sa.WSTable(raster_weight_list)
        self.suitability_map = arcpy.sa.WeightedSum(WSumTableObj)

        x = []
        for r, c in self.suitability_map:
            x.append(self.suitability_map[r, c])
        self.suitability_values = list(filter(lambda v: not math.isnan(v), x))

    def show_stats(self):
        print('Mean: {}'.format(self.suitability_map.mean))
        print('Min: {}'.format(self.suitability_map.minimum))
        print('Max: {}'.format(self.suitability_map.maximum))

    def show_hist(self, n_bins=20):
        fig, ax1 = plt.subplots()
        ax1.hist(self.suitability_values, bins=n_bins)
        ax1.set_xlabel('Suitability Map')
        ax1.set_ylabel('Count')
        plt.title('Histogram of Suitability Map')
        plt.show()


def main():
    # test code here
    c1 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\dem_24'))
    c2 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\dem_24'))
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
    c1.transform('range', transform_params_range)
    transform_params_continous = {
        'name': 'small',
        'from_scale': 1,
        'to_scale': 10
    }
    c2.transform('continous', transform_params_continous)
    s = SuitabilityModel()
    s.add_criteria(c1)
    s.add_criteria(c2)
    s.calculate()
    s.show_stats()
    s.show_hist()


if __name__ == "__main__":
    main()
