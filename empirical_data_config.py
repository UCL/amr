# Empirical Data Source Configuration
# Update this file with actual API endpoints, credentials, and data access details

import os
from typing import Dict, Any

class DataSourceConfig:
    """Configuration for empirical data sources"""
    
    # ECDC (European Centre for Disease Prevention and Control)
    ECDC_CONFIG = {
        'base_url': 'https://www.ecdc.europa.eu/en/antimicrobial-resistance/surveillance-and-disease-data/data-ecdc',
        'api_key': os.getenv('ECDC_API_KEY', ''),  # Set environment variable
        'datasets': {
            'resistance': 'AMR_surveillance_data',
            'consumption': 'antimicrobial_consumption_ESAC-Net'
        },
        'file_paths': {
            # Download these files manually and update paths
            'resistance_csv': './data/ecdc_resistance_2023.csv',
            'consumption_csv': './data/ecdc_consumption_2023.csv'
        }
    }
    
    # WHO GLASS (Global Antimicrobial Resistance Surveillance System)
    WHO_GLASS_CONFIG = {
        'base_url': 'https://www.who.int/data/gho/data/themes/world-health-statistics',
        'api_key': os.getenv('WHO_API_KEY', ''),
        'datasets': {
            'resistance_rates': 'GLASS_AMR_data',
            'surveillance_coverage': 'GLASS_country_participation'
        },
        'file_paths': {
            # WHO provides Excel/CSV downloads
            'glass_data_2022': './data/who_glass_2022.xlsx',
            'glass_data_2021': './data/who_glass_2021.xlsx'
        }
    }
    
    # IQVIA (Commercial pharmaceutical data)
    IQVIA_CONFIG = {
        'base_url': 'https://www.iqvia.com/locations/united-states/solutions/commercial-effectiveness/brand-strategy-and-management/midas',
        'api_key': os.getenv('IQVIA_API_KEY', ''),  # Requires commercial license
        'credentials': {
            'username': os.getenv('IQVIA_USERNAME', ''),
            'password': os.getenv('IQVIA_PASSWORD', '')
        },
        'datasets': {
            'midas_sales': 'pharmaceutical_sales_by_country',
            'market_share': 'antibiotic_market_penetration'
        },
        'file_paths': {
            # IQVIA provides custom data exports
            'sales_data_2023': './data/iqvia_sales_2023.csv',
            'market_data_2023': './data/iqvia_market_2023.csv'
        }
    }
    
    # CDC (Centers for Disease Control and Prevention)
    CDC_CONFIG = {
        'base_url': 'https://www.cdc.gov/drugresistance/biggest-threats.html',
        'api_endpoints': {
            'nndss': 'https://wwwn.cdc.gov/nndss/data-and-statistics.html',
            'ar_threats': 'https://www.cdc.gov/drugresistance/biggest-threats.html'
        },
        'file_paths': {
            # CDC provides public data downloads
            'ar_threats_2019': './data/cdc_ar_threats_2019.pdf',  # May need manual extraction
            'nndss_data': './data/cdc_nndss_2022.csv'
        }
    }
    
    # National Statistics Offices
    NATIONAL_STATS_CONFIG = {
        'sources': {
            'usa': {
                'url': 'https://www.cdc.gov/nchs/nvss/deaths.htm',
                'mortality_data': './data/usa_mortality_icd10.csv'
            },
            'uk': {
                'url': 'https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths',
                'mortality_data': './data/uk_mortality_2022.csv'
            },
            'germany': {
                'url': 'https://www.destatis.de/EN/Themes/Society-Environment/Population/Deaths-Life-Expectancy/_node.html',
                'mortality_data': './data/germany_mortality_2022.csv'
            }
        }
    }
    
    # Regional mappings for data harmonization
    REGION_MAPPINGS = {
        'countries_to_regions': {
            # North America
            'united_states': 'north_america',
            'canada': 'north_america',
            'mexico': 'north_america',
            
            # Europe
            'germany': 'europe',
            'france': 'europe',
            'united_kingdom': 'europe',
            'italy': 'europe',
            'spain': 'europe',
            'netherlands': 'europe',
            'sweden': 'europe',
            'norway': 'europe',
            
            # Asia
            'china': 'asia',
            'india': 'asia',
            'japan': 'asia',
            'south_korea': 'asia',
            'thailand': 'asia',
            'indonesia': 'asia',
            
            # Africa
            'south_africa': 'africa',
            'nigeria': 'africa',
            'kenya': 'africa',
            'egypt': 'africa',
            
            # South America
            'brazil': 'south_america',
            'argentina': 'south_america',
            'colombia': 'south_america',
            'chile': 'south_america',
            
            # Oceania
            'australia': 'oceania',
            'new_zealand': 'oceania'
        },
        
        'ecdc_countries': [
            'austria', 'belgium', 'bulgaria', 'croatia', 'cyprus', 'czech_republic',
            'denmark', 'estonia', 'finland', 'france', 'germany', 'greece',
            'hungary', 'iceland', 'ireland', 'italy', 'latvia', 'liechtenstein',
            'lithuania', 'luxembourg', 'malta', 'netherlands', 'norway', 'poland',
            'portugal', 'romania', 'slovakia', 'slovenia', 'spain', 'sweden', 'united_kingdom'
        ]
    }
    
    # Data quality indicators
    QUALITY_THRESHOLDS = {
        'min_sample_size': 100,        # Minimum specimens for resistance data
        'min_country_coverage': 0.7,   # Minimum proportion of countries reporting
        'max_missing_years': 3,        # Maximum consecutive missing years
        'confidence_level': 0.95       # Statistical confidence level
    }

# Environment setup instructions
SETUP_INSTRUCTIONS = """
EMPIRICAL DATA SOURCE SETUP INSTRUCTIONS
=======================================

1. ECDC Data:
   - Visit: https://www.ecdc.europa.eu/en/antimicrobial-resistance/surveillance-and-disease-data/data-ecdc
   - Download: "Antimicrobial resistance surveillance in Europe" annual reports
   - Save CSV files to ./data/ directory

2. WHO GLASS Data:
   - Visit: https://www.who.int/initiatives/glass
   - Download: GLASS surveillance reports (Excel format)
   - Convert to CSV and save to ./data/ directory

3. IQVIA Data (Commercial License Required):
   - Contact IQVIA for MIDAS database access
   - Set environment variables: IQVIA_API_KEY, IQVIA_USERNAME, IQVIA_PASSWORD
   - Configure data export paths

4. CDC Data:
   - Visit: https://www.cdc.gov/drugresistance/biggest-threats.html
   - Download: AR Threats Report data
   - Visit: https://wwwn.cdc.gov/nndss/ for surveillance data

5. National Statistics:
   - Download mortality data from national statistics offices
   - Convert to standardized CSV format

6. Environment Variables:
   Set these in your environment:
   export ECDC_API_KEY="your_key_here"
   export WHO_API_KEY="your_key_here"  
   export IQVIA_API_KEY="your_key_here"
   export IQVIA_USERNAME="your_username"
   export IQVIA_PASSWORD="your_password"

7. Data Directory Structure:
   Create: ./data/
   ├── ecdc_resistance_2023.csv
   ├── ecdc_consumption_2023.csv
   ├── who_glass_2022.xlsx
   ├── iqvia_sales_2023.csv
   ├── cdc_ar_threats_2019.csv
   └── national_mortality_data/
       ├── usa_mortality_icd10.csv
       ├── uk_mortality_2022.csv
       └── germany_mortality_2022.csv
"""
