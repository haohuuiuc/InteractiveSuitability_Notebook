import arcpy
from Criteria import *


class SuitabilityModel:
    def __init__(self, weight_method='multiplier', from_scale=1, to_scale=10):
        self.criteria = []
        self.weight = []
        self.from_scale = from_scale
        self.to_scale = to_scale
        self.suitability_map = None

    def add_criteria(self, criterion, weight=1):
        self.criteria.append(criterion)
        self.weight.append(weight)

    def calculate(self):
        [["snow", "VALUE", 0.25], ["land", "VALUE", 0.25],
         ["soil", "VALUE", 0.5]]
        raster_weight_list = []
        for i in range(len(self.criteria)):
            raster_weight_list.append([self.criteria[i].raster, 'VALUE', self.weight[i]])
        WSumTableObj  = arcpy.sa.WSTable(raster_weight_list)
        self.suitability_map = arcpy.sa.WeightedSum(WSumTableObj)

    def stats(self):
        x = []
        for r, c in self.suitability_map:
            x.append(self.suitability_map[r, c])
        self.suitability_values = list(filter(lambda v: not math.isnan(v), x))
        print('Mean: {}'.format(statistics.mean(self.suitability_values)))
        print('Min: {}'.format(min(self.suitability_values)))
        print('Max: {}'.format(max(self.suitability_values)))


def main():
    # test code here
    c1 = Criteria(arcpy.Raster(r'\\archive\CRData\ArcGISPro\raster-analysis\SuitabilityData\dem_24'))
    s = SuitabilityModel()
    s.add_criteria(c1)
    s.calculate()
    s.stats()

if __name__ == "__main__":
    main()
