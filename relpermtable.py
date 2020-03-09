import pandas as pd
import numpy as np
import time
# =============================================================================
# Class constructor that reads the ASCII file and calculates unscaled
# saturation endpoints
# =============================================================================
class RelPermTable:
    def __init__(self, fol_path):
        self.df = pd.DataFrame(columns=['Sg', 'Krg', 'Krw', 'Pc']) # dataframe storing original unscaled values
#        self.sc_df = pd.DataFrame(columns=['Sg', 'Sw', 'temp_S', 'Krg', 'Krw', 'Pc'])
#       Line counter to know how many lines are in the original relative permeability table
        self.line_counter = 0
#       Open the file and populate the dataframe with relative permeability points
        f = open(fol_path, 'r')
        for line in f:
            if 'SGWFN' in line:
                for line in f:
                    line = line.strip()
                    if (line.split() and not line.startswith('/') and not line.startswith('--')):
                        cols = line.split()
                        self.df = self.df.append(pd.Series(([float(i) for i in cols]), index=['Sg', 'Krg', 'Krw', 'Pc']), ignore_index=True)
                        self.line_counter += 1
#       Initialise unscaled saturation end points
        self.sgmax = self.df['Sg'].max()
        self.sgco = self.df['Sg'].min()

        self.swco = 1 - self.sgmax
        self.swmax = 1- self.swco
    
        self.sgcr = self.df.loc[self.df['Krg'].eq(0), 'Sg'].tail(1).iloc[0]
        self.swcr = 1 - self.df.loc[self.df['Krw'].eq(0), 'Sg'].head(1).iloc[0]
        
        self.sg_krgr = 1 - self.swcr
        self.sw_krwr = 1 - self.sgcr
        
        self.krg_max = self.df.Krg.max()
        self.krw_max = self.df.Krw.max()
        
        self.krwr = self.df.loc[np.isclose(1 - self.df.Sg, self.sw_krwr), 'Krw'].iloc[0]
        self.krgr = self.df.loc[np.isclose(self.df.Sg, self.sg_krgr), 'Krg'].iloc[0]
        
    def do_save_data(self, fol_path):
        self.sc_df.to_csv(fol_path, sep = '\t', float_format='%.8f', columns = ['Sg', 'Sw', 'Krg', 'Krw', 'Pc'])
        
# =============================================================================
#     Method to insert saturation endpoints into fishbone saturation array
# =============================================================================
    def __insert(self, array, value):
        array = np.asarray(array)
        if not np.any(np.isclose(array, value)):
            array = np.insert(array, np.searchsorted(array, value), value)
        return array

# =============================================================================
#     Method to do horizontal and vertical saturation and relative
#     permeability scaling
# =============================================================================
    def three_point_scaling(self, set_SWCR = None, set_SWU = None, 
                            set_KRWR = None, set_KRW = None, 
                            set_SGCR = None, set_SGU = None,  
                            set_KRGR = None, set_KRG = None):
#        print('Class: ', set_SWU, set_SGU, set_SGCR, set_SWCR)
        if set_SWCR is None: set_SWCR = self.swcr
        if set_SWU is None: set_SWU = self.swmax
        if set_KRWR is None: set_KRWR = self.krwr
        if set_KRW is None: set_KRW = self.krw_max
        if set_SGCR is None: set_SGCR = self.sgcr
        if set_SGU is None: set_SGU = self.sgmax
        if set_KRGR is None: set_KRGR = self.krgr
        if set_KRGR is None: set_KRG = self.krg_max
        self.sc_df = pd.DataFrame(columns=['Sg', 'Sw', 'temp_S', 'Krg', 'Krw', 'Pc'])
        set_SGL = 1 - set_SWU
        set_SG_KRGR = 1 - set_SWCR
        set_SW_KRWR = 1 - set_SGCR
        
#       Generating the fishbone for scaled relative permeability
#       equally spaced between 0 and 1
        arr = pd.Series(np.linspace(0,1,21)) 
        
        # Adding SWL, SGL, SWCR and SGCR
        arr = self.__insert(arr, set_SGL)
        arr = self.__insert(arr, set_SGU)
        arr = self.__insert(arr, set_SGCR)
        arr = self.__insert(arr, set_SG_KRGR)
        self.sc_df['Sg'] = arr
        self.sc_df['Sw'] = 1 - self.sc_df['Sg']

#       Filling scaled Krg values
        self.sc_df.loc[self.sc_df.Sg <= set_SGCR, 'Krg'] = 0
        self.sc_df.loc[self.sc_df.Sg >= set_SGU, 'Krg'] = self.krg_max

        y = lambda x: self.sgcr + ((x - set_SGCR)*(self.sg_krgr - self.sgcr)/
                                   (set_SG_KRGR - set_SGCR)) if x <= set_SG_KRGR else self.sg_krgr + ((x - set_SG_KRGR)*
                                   (self.sgmax - self.sg_krgr)/
                                   (set_SGU-set_SG_KRGR))

#        temp_S is a pseudo gas saturation column used for looking up scaled Krg
        self.sc_df['temp_S'] = self.sc_df['Sg'].apply(y)
        
#        Lookup Krg based on pseudo gas saturation Sg'
        self.sc_df.loc[(self.sc_df.Sg < set_SGU) & (self.sc_df.Sg > set_SGCR), 'Krg'] = np.interp(
                self.sc_df.loc[(self.sc_df.Sg < set_SGU) & (self.sc_df.Sg > set_SGCR), 'temp_S'].values, 
                                     self.df['Sg'].values, self.df['Krg'].values)

#       Filling scaled Krw values
        self.sc_df.loc[self.sc_df.Sw <= set_SWCR, 'Krw'] = 0
        self.sc_df.loc[self.sc_df.Sw >= set_SWU, 'Krw'] = self.krw_max
        
        yw = lambda xw: self.swcr + ((xw - set_SWCR)*(self.sw_krwr - self.swcr)/
                                   (set_SW_KRWR - set_SWCR)) if xw <= set_SW_KRWR else self.sw_krwr + ((xw - set_SW_KRWR)*
                                   (self.swmax - self.sw_krwr)/(set_SWU-set_SW_KRWR))
#       Water look up saturation values are converted to gas saturation
#       to be able to interpolate through the original relative permeability
#       table and to not calculate there water saturation in the original table
        self.sc_df['temp_S'] = 1-self.sc_df['Sw'].apply(yw)
#       NOTE: !!! Interpolation look up is done on gas saturation - it is correct as temp_S is 
#       converted from water saturation from gas. All done to not add an additional
#       column to the original relative permeability dataframe
        self.sc_df.loc[(self.sc_df.Sw < set_SWU) & (self.sc_df.Sw > set_SWCR), 'Krw'] = np.interp(
                self.sc_df.loc[(self.sc_df.Sw < set_SWU) & (self.sc_df.Sw > set_SWCR), 'temp_S'].values, 
                                     self.df['Sg'].values, self.df['Krw'].values)
        
#       Vertical scaling part of the method
        wat = lambda x: x['Krw'] * set_KRWR / self.krwr if ((x['Sw'] <= set_SW_KRWR) & 
                             (x['Sw'] > set_SWCR)) else (0 if x['Sw'] <= set_SWCR else 
                             ((set_KRWR + ((x['Krw'] - self.krwr)*(set_KRW - set_KRWR)/(self.krw_max - self.krwr))) if self.krw_max != self.krwr else np.NaN) )
        self.sc_df['Krw'] = self.sc_df.apply(wat, axis = 1)
        
        if self.sc_df['Krw'].isnull().values.any():
            self.sc_df.loc[np.isclose(self.sc_df.Sw, set_SW_KRWR), 'Krw'] = set_KRWR
            self.sc_df.loc[np.isclose(self.sc_df.Sw, self.sc_df.Sw.max()), 'Krw'] = set_KRW
            self.sc_df['Krw'] = self.sc_df.interpolate(method = 'linear', limit_area = 'inside')['Krw']
        
        gas = lambda x: x['Krg'] * set_KRGR / self.krgr if ((x['Sg'] <= set_SG_KRGR) & 
                             (x['Sg'] > set_SGCR)) else (0 if x['Sg'] <= set_SGCR else 
                             ((set_KRGR + ((x['Krg'] - self.krgr)*(set_KRG - set_KRGR)/(self.krg_max - self.krgr))) if self.krg_max != self.krgr else np.NaN) )
        self.sc_df['Krg'] = self.sc_df.apply(gas, axis = 1)
        

        if self.sc_df['Krg'].isnull().values.any():
            self.sc_df.loc[np.isclose(self.sc_df.Sg, set_SG_KRGR), 'Krg'] = set_KRGR
            self.sc_df.loc[np.isclose(self.sc_df.Sg, self.sc_df.Sg.max()), 'Krg'] = set_KRG
            self.sc_df['Krg'] = self.sc_df.interpolate(method = 'linear', limit_area = 'inside')['Krg']
            
        self.sc_df['temp_S'] = np.nan
        pc = lambda x: self.sgco + ((x - set_SGL)*(self.sgmax - self.sgco)/(set_SGU - set_SGL))
        self.sc_df['temp_S'] = self.sc_df['Sg'].apply(pc)
        self.sc_df['Pc'] = np.interp(
                self.sc_df['temp_S'].values, 
                                     self.df['Sg'].values, self.df['Pc'].values)
            
#        self.sc_df = self.sc_df.drop('temp_S', axis=1) 
#        self.sc_df[self.sc_df['Sg'] >= set_SGL]

if __name__ == '__main__':
    start_time = time.time()
    object = RelPermTable(".\DATA\EPS.INC")

    set_SWCR = 0.15
    set_SGCR = 0.36
    set_SWU = 0.95
    set_SGU = 1.0
#    
    set_KRG = 1.0
    set_KRW = 0.65
    set_KRGR = 0.9
    set_KRWR = 0.3
#    for i in np.linspace(0.3, 0.4, 3):
#   set_SWCR, set_SWU, set_KRWR, set_KRW, set_SGCR, set_SGU,  set_KRGR, set_KRG
    object.three_point_scaling(set_SWCR, set_SWU, set_KRWR, set_KRW, set_SGCR, set_SGU, set_KRGR, set_KRG)
#    object.three_point_scaling()
#    
    print(object.sc_df)
    print('\nReduced:')
#    print(object.sc_df[np.isclose(object.sc_df['Sg']-(1 - set_SWU), 0) or (object.sc_df['Sg'] > (1 - set_SWU))])
    print(object.sc_df[object.sc_df['Sg'].values > 0.05])
#    print(object.sc_df[abs(object.sc_df['Sg'] - (1 - set_SWU)) < 1e-15])
#    print(object.sc_df[~(object.sc_df['Sg'].isin(object.sc_df[object.sc_df['Sg'].values < (1 - set_SWU)]['Sg']))])
    print('\nUnscaled Swco: ', object.swco)
    print('Unscaled Swcr: ', object.swcr)
    print('Unscaled Sw(krwr): ', object.sw_krwr)
    print('Unscaled Sw max: ', object.swmax)
    print('\nUnscaled krwr: ', object.krwr)
    print('Unscaled krw max: ', object.krw_max)
    print('\nUnscaled Sgco: ', object.sgco)
    print('Unscaled Sgcr: ', object.sgcr)
    print('Unscaled Sg(krgr): ', object.sg_krgr)
    print('Unscaled Sg max: ', object.sgmax)
    print('\nUnscaled krgr: ', object.krgr)
    print('Unscaled krg max: ', object.krg_max)
    print('\nIt took', time.time()-start_time, 'seconds.')   
    
    
#    print('Set SWCR: ', set_SWCR)
#    print('Set SGCR: ', set_SGCR)
#    print('Set SWU: ', set_SWU)
    
#    set_SWCR, set_KRWR, set_KRW, , set_KRGR, set_KRG, set_SGCR
#    object.three_point_vert_scaling(set_SWCR, 0.3, new_KRW, 0.9, new_KRG, set_SGCR)
  

    
#    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#        print(new_df[['sc_Sg', 'ps_Sg']])
#    print(" swco =", new_swco, "\n", "sgmax =", new_sgmax, "\n", "sgco =", 
#          new_sgco, "\n", "swmax =", new_swmax,
#          "\n", "sgcr =", new_sgcr, "\n", "swcr =", new_swcr)
#    for index, row in new_df.iterrows():
#        print(type(index), type(row["Sg"]))