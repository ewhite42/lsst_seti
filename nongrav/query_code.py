## This program is meant to retrieve the non-gravitational 
## acceleration parameters of several comets stored in the 
## CODE catalog: https://pad2.astro.amu.edu.pl/comets/index.php

## MEW, 23 Aug 2025

import pylab as plt
import pandas as pd
import requests
from bs4 import BeautifulSoup as bs

def main():
    ## first, read in a list of the objects with nongrav
    ## acceleration parameters:
    df_ng = pd.read_csv('/home/ellie/Downloads/codec_results.txt', sep='\\s+', engine='python')
    df_ng.columns = ['designation', 'model', 'nobs', 'arc1', 'arc2', 'epoch', \
                     'T[TT]', 'q[au]', 'e', 'omega [deg]', 'ohm [deg]', \
                     'i[deg]', '1/a', 'dT', 'dq', 'de', 'domega', \
                     'dOmega', 'di', 'drecip(a)']
                          
    #print(df_ng.head(1))
    
    ## next, query the page that lists rows for all available
    ## comets, then find which ones match those in df_ng
    url = 'https://pad2.astro.amu.edu.pl/comets/index.php'
    df_idx = pd.read_html(url)
    df_idx = df_idx[0]
    
    df_idx.columns = ['designation', 'model', 'name', 'class', 'observations',\
                      'date interval', 'epoch', 'T[TT]', 'q[au]', 'e', 'omega [deg]',\
                      'ohm [deg]', 'i[deg]', '1/a']
                      
    #print(df_idx.head(1))
    #print(df_idx.columns)
    
    response = requests.get(url)
    f = response.text
    
    soup = bs(f, 'lxml')
    parsed_table = soup.find_all('table')[0] 
    data = [[td.a['href'] if td.find('a') else 
             ''.join(td.stripped_strings)
             for td in row.find_all('td')]
            for row in parsed_table.find_all('tr')]
    df_html = pd.DataFrame(data[1:])
    df_html.columns = ['url', 'model', 'name', 'class', 'observations',\
                      'date interval', 'epoch', 'T[TT]', 'q[au]', 'e', 'omega [deg]',\
                      'ohm [deg]', 'i[deg]', '1/a']
    #print(df_html)
    
    designations = df_idx['designation']
    df_html.drop(index=0, inplace=True)
    df_html.reset_index(drop=True, inplace=True)
    df_html.insert(loc=0, column='designation', value=designations)
    #print(df_html.head(1))   
    
    url_prefix = 'https://pad2.astro.amu.edu.pl/comets/'
    nongrav_rows = []
    counter = 0
    
    for index, row in df_html.iterrows():
        comet_url = url_prefix+str(row['url'])
        
        ## read in the tables from the page for each comet
        dfs = pd.read_html(comet_url)
        #print(len(dfs))
        
        if len(dfs) == 3:
            #print(tables[2])
        
            ## make the table of non-gravitational parameters
            df = dfs[-1].transpose()
            df.reset_index(drop=True, inplace=True)
            df.columns = df.iloc[0]
            df.drop(index=0, inplace=True)
            df.reset_index(drop=True, inplace=True)
                
            df_html.at[index, 'nongrav bool'] = True
            df_html.at[index, 'A1 [10-8au/day2]'] = df.at[0, 'A1 [10-8au/day2]']
            df_html.at[index, 'A2 [10-8au/day2]'] = df.at[0, 'A2 [10-8au/day2]']
            df_html.at[index, 'A3 [10-8au/day2]'] = df.at[0, 'A3 [10-8au/day2]']
            
            counter += 1
            print('{} % complete'.format(round(100*counter/len(df_html), 3)))
            
        else:
            df_html.at[index, 'nongrav bool'] = False
            counter += 1
            print('{} % complete'.format(round(100*counter/len(df_html), 3)))
            continue 
        
    print(len(df_html))
    df_nongrav = df_html[df_html['nongrav bool'] == True]
    print(len(df_nongrav))
    
    nbins = 20
    plt.hist(df_nongrav['A1 [10-8au/day2]'], nbins)
    plt.show()
    
    plt.hist(df_nongrav['A2 [10-8au/day2]'], nbins)
    plt.show()
    
    plt.hist(df_nongrav['A3 [10-8au/day2]'], nbins)
    plt.show()
    
if __name__ == '__main__':
    main()
    
