# Enhanced Empirical Data Configuration
# Comprehensive configuration for acquiring real consumption data across all regions and drugs

import os
from typing import Dict, Any, List

class EnhancedDataSourceConfig:
    """Enhanced configuration for comprehensive empirical data acquisition"""
    
    # === PRIORITY DATA SOURCES (Immediate Implementation) ===
    
    # ECDC - European Centre for Disease Prevention and Control
    ECDC_ENHANCED_CONFIG = {
        'base_url': 'https://www.ecdc.europa.eu/en/antimicrobial-consumption/database',
        'esac_net_url': 'https://www.ecdc.europa.eu/sites/default/files/documents/ESAC-Net_Report_2022.xlsx',
        'file_paths': {
            'consumption_2022': './data/ecdc_consumption_2022_complete.xlsx',
            'consumption_2021': './data/ecdc_consumption_2021_complete.xlsx', 
            'consumption_2020': './data/ecdc_consumption_2020_complete.xlsx',
            'quality_indicators': './data/ecdc_quality_indicators_2022.xlsx'
        },
        'drugs_covered': [
            # Beta-lactams
            'penicilling',         # J01CA (Penicillins with extended spectrum)
            'ampicillin',          # J01CA01
            'amoxicillin',         # J01CA04
            'amoxicillin_clavulanate',  # J01CR02
            'piperacillin_tazobactam',  # J01CR05
            'cephalexin',          # J01DB01
            'cefuroxime',          # J01DC02
            'ceftriaxone',         # J01DD04
            'ceftazidime',         # J01DD02
            'cefepime',            # J01DE01
            'meropenem',           # J01DH02
            'imipenem_c',          # J01DH51
            'ertapenem',           # J01DH03
            
            # Macrolides
            'erythromycin',        # J01FA01
            'azithromycin',        # J01FA10
            'clarithromycin',      # J01FA09
            'clindamycin',         # J01FF01
            
            # Fluoroquinolones
            'ciprofloxacin',       # J01MA02
            'levofloxacin',        # J01MA12
            'moxifloxacin',        # J01MA14
            'ofloxacin',           # J01MA01
            
            # Tetracyclines
            'tetracycline',        # J01AA02
            'doxyclycline',        # J01AA02
            
            # Aminoglycosides
            'gentamicin',          # J01GB03
            'tobramycin',          # J01GB01
            'amikacin',            # J01GB06
            
            # Others
            'vancomycin',          # J01XA01
            'linezolid',           # J01XX08
            'trim_sulf',           # J01EE01
            'nitrofurantoin',      # J01XE01
            'metronidazole',       # J01XD01
            'rifampicin'           # J04AB02
        ],
        'countries_covered': [
            'austria', 'belgium', 'bulgaria', 'croatia', 'cyprus', 'czech_republic',
            'denmark', 'estonia', 'finland', 'france', 'germany', 'greece',
            'hungary', 'iceland', 'ireland', 'italy', 'latvia', 'liechtenstein',
            'lithuania', 'luxembourg', 'malta', 'netherlands', 'norway', 'poland',
            'portugal', 'romania', 'slovakia', 'slovenia', 'spain', 'sweden', 'united_kingdom'
        ],
        'data_format': 'DDD per 1000 inhabitants per day',
        'conversion_factor': 36.5,  # Convert to courses per 100k per year
        'quality_threshold': 80  # Minimum data quality percentage
    }
    
    # WHO Global Health Observatory
    WHO_GHO_CONFIG = {
        'base_url': 'https://www.who.int/data/gho/data/themes/world-health-statistics',
        'api_url': 'https://ghoapi.azureedge.net/api/',
        'endpoints': {
            'antimicrobial_consumption': 'WHS4_543',  # AMC indicator
            'antimicrobial_resistance': 'AMR_GLASS',   # GLASS surveillance
        },
        'file_paths': {
            'gho_consumption_global': './data/who_gho_consumption_global.json',
            'gho_consumption_regional': './data/who_gho_consumption_regional.json'
        },
        'regions_mapping': {
            'AFR': 'africa',
            'AMR': 'north_america',
            'SEAR': 'asia', 
            'EUR': 'europe',
            'EMR': 'asia',  # Eastern Mediterranean -> Asia
            'WPR': 'asia'   # Western Pacific -> Asia
        }
    }
    
    # Australian AURA Surveillance System
    AUSTRALIA_AURA_CONFIG = {
        'base_url': 'https://www.safetyandquality.gov.au/our-work/antimicrobial-resistance',
        'data_url': 'https://www.safetyandquality.gov.au/sites/default/files/2023-06/aura-2022-report.pdf',
        'file_paths': {
            'aura_2022': './data/australia_aura_2022.pdf',
            'aura_2021': './data/australia_aura_2021.pdf',
            'aura_consumption_data': './data/australia_consumption_extracted.csv'
        },
        'drugs_covered': [
            'penicilling', 'amoxicillin', 'amoxicillin_clavulanate',
            'cephalexin', 'cefuroxime', 'ceftriaxone', 'cefepime',
            'erythromycin', 'azithromycin', 'clarithromycin', 'clindamycin',
            'ciprofloxacin', 'levofloxacin', 'tetracycline', 'doxyclycline',
            'gentamicin', 'vancomycin', 'linezolid', 'trim_sulf',
            'nitrofurantoin', 'metronidazole'
        ],
        'data_extraction_required': True,  # PDF extraction needed
        'extraction_tool': 'tabula-py'
    }
    
    # === COMMERCIAL DATA SOURCES ===
    
    # Enhanced IQVIA Configuration
    IQVIA_ENHANCED_CONFIG = {
        'products': {
            'midas': {
                'description': 'Global pharmaceutical sales database',
                'coverage': 'Global',
                'cost': '$75,000-200,000 annually',
                'contact': 'https://www.iqvia.com/solutions/commercialization/brand-strategy-and-management/midas'
            },
            'national_sales_perspectives': {
                'description': 'Country-specific pharmaceutical sales',
                'regions': ['USA', 'Canada', 'Germany', 'France', 'UK', 'Japan', 'Australia'],
                'cost': '$25,000-50,000 per country annually'
            }
        },
        'alternative_vendors': {
            'pharmaintelligence': {
                'url': 'https://www.pharmaintelligence.informa.com/',
                'products': ['Datamonitor Healthcare', 'Global Data']
            },
            'evaluate_pharma': {
                'url': 'https://www.evaluate.com/',
                'focus': 'Pharmaceutical market intelligence'
            }
        }
    }
    
    # === NATIONAL DATA SOURCES ===
    
    NATIONAL_SOURCES_CONFIG = {
        'asia': {
            'japan': {
                'agency': 'National Institute of Infectious Diseases (NIID)',
                'url': 'https://www.niid.go.jp/niid/en/',
                'data_source': 'Japan Nosocomial Infections Surveillance (JANIS)',
                'file_path': './data/japan_niid_consumption.xlsx',
                'contact': 'janis@niid.go.jp'
            },
            'south_korea': {
                'agency': 'Korea Disease Control and Prevention Agency (KDCA)',
                'url': 'https://www.kdca.go.kr/',
                'data_source': 'Korean Antimicrobial Resistance Monitoring System',
                'file_path': './data/korea_karms_consumption.xlsx'
            },
            'china': {
                'agency': 'National Health Commission',
                'note': 'Limited public access - requires academic collaboration',
                'potential_contact': 'Chinese Center for Disease Control and Prevention'
            },
            'india': {
                'agency': 'Indian Council of Medical Research (ICMR)',
                'url': 'https://www.icmr.gov.in/',
                'data_source': 'AMR surveillance network',
                'file_path': './data/india_icmr_consumption.pdf'
            }
        },
        
        'africa': {
            'south_africa': {
                'agency': 'National Institute for Communicable Diseases (NICD)',
                'url': 'https://www.nicd.ac.za/',
                'data_source': 'SASCM surveillance',
                'file_path': './data/south_africa_nicd_consumption.pdf'
            },
            'nigeria': {
                'agency': 'Nigeria Centre for Disease Control (NCDC)',
                'url': 'https://ncdc.gov.ng/',
                'file_path': './data/nigeria_ncdc_amr_report.pdf'
            },
            'kenya': {
                'agency': 'Kenya Medical Research Institute (KEMRI)',
                'url': 'https://www.kemri.org/',
                'file_path': './data/kenya_kemri_amr_data.xlsx'
            }
        },
        
        'south_america': {
            'brazil': {
                'agency': 'National Health Surveillance Agency (ANVISA)',
                'url': 'https://www.gov.br/anvisa/',
                'data_source': 'Antibiotic consumption monitoring',
                'file_path': './data/brazil_anvisa_consumption.xlsx'
            },
            'argentina': {
                'agency': 'National Administration of Medicines (ANMAT)',
                'url': 'https://www.argentina.gob.ar/anmat',
                'file_path': './data/argentina_anmat_antibiotics.pdf'
            },
            'colombia': {
                'agency': 'National Institute for Drug and Food Surveillance (INVIMA)',
                'url': 'https://www.invima.gov.co/',
                'file_path': './data/colombia_invima_consumption.xlsx'
            }
        }
    }
    
    # === IMPLEMENTATION PRIORITIES ===
    
    IMPLEMENTATION_PRIORITY = [
        {
            'priority': 1,
            'source': 'ECDC_ENHANCED',
            'expected_improvement': '+22 drugs for Europe',
            'effort': 'Low - direct download',
            'timeline': '1-2 weeks'
        },
        {
            'priority': 2,
            'source': 'AUSTRALIA_AURA',
            'expected_improvement': '+20 drugs for Oceania',
            'effort': 'Medium - PDF extraction',
            'timeline': '2-3 weeks'
        },
        {
            'priority': 3,
            'source': 'WHO_GHO_API',
            'expected_improvement': '+8 drugs for 3 regions',
            'effort': 'Medium - API integration',
            'timeline': '3-4 weeks'
        },
        {
            'priority': 4,
            'source': 'NATIONAL_SOURCES',
            'expected_improvement': '+15 drugs for Asia/Africa/S.America',
            'effort': 'High - manual extraction',
            'timeline': '6-8 weeks'
        },
        {
            'priority': 5,
            'source': 'IQVIA_EXPANSION',
            'expected_improvement': '+35 drugs for all regions',
            'effort': 'High - commercial negotiation',
            'timeline': '3-6 months'
        }
    ]
    
    # === DRUG MAPPING FOR HARMONIZATION ===
    
    DRUG_HARMONIZATION_MAP = {
        # ATC Code -> Simulation Drug Name
        'J01CA04': 'amoxicillin',
        'J01CR02': 'amoxicillin_clavulanate',
        'J01CA01': 'ampicillin',
        'J01CR01': 'ampicillin_sulbactam',
        'J01CR05': 'piperacillin_tazobactam',
        'J01DB01': 'cephalexin',
        'J01DC02': 'cefuroxime',
        'J01DD04': 'ceftriaxone',
        'J01DD02': 'ceftazidime',
        'J01DE01': 'cefepime',
        'J01DH02': 'meropenem',
        'J01DH51': 'imipenem_c',
        'J01DH03': 'ertapenem',
        'J01FA01': 'erythromycin',
        'J01FA10': 'azithromycin',
        'J01FA09': 'clarithromycin',
        'J01FF01': 'clindamycin',
        'J01MA02': 'ciprofloxacin',
        'J01MA12': 'levofloxacin',
        'J01MA14': 'moxifloxacin',
        'J01MA01': 'ofloxacin',
        'J01AA02': 'tetracycline',
        'J01AA12': 'doxyclycline',
        'J01GB03': 'gentamicin',
        'J01GB01': 'tobramycin',
        'J01GB06': 'amikacin',
        'J01XA01': 'vancomycin',
        'J01XX08': 'linezolid',
        'J01EE01': 'trim_sulf',
        'J01XE01': 'nitrofurantoin',
        'J01XD01': 'metronidazole',
        'J04AB02': 'rifampicin'
    }

# Setup instructions for enhanced data acquisition
ENHANCED_SETUP_INSTRUCTIONS = """
ENHANCED EMPIRICAL DATA ACQUISITION SETUP
========================================

PHASE 1: FREE/PUBLIC SOURCES (1-4 weeks)
----------------------------------------

1. ECDC ESAC-Net Data:
   a) Visit: https://www.ecdc.europa.eu/en/antimicrobial-consumption/database
   b) Download complete surveillance reports (2020-2023)
   c) Files needed:
      - ESAC-Net_Report_2022.xlsx (latest complete data)
      - Country-specific consumption tables
      - Data quality indicators

2. WHO Global Health Observatory:
   a) API access: https://ghoapi.azureedge.net/api/
   b) Register for API key (free)
   c) Query indicators: WHS4_543 (consumption), AMR_GLASS (resistance)

3. Australia AURA Surveillance:
   a) Download: https://www.safetyandquality.gov.au/sites/default/files/2023-06/aura-2022-report.pdf
   b) Install tabula-py for PDF extraction: pip install tabula-py
   c) Extract consumption tables from report

4. National Surveillance Reports:
   - Japan NIID: Download annual AMR surveillance reports
   - South Korea KDCA: Access KARMS consumption data
   - Brazil ANVISA: Antibiotic consumption monitoring reports

PHASE 2: ACADEMIC PARTNERSHIPS (2-6 months)
------------------------------------------

1. Contact WHO AMR Secretariat:
   Email: amr@who.int
   Request: Access to unpublished GLASS consumption data

2. Academic collaborations:
   - London School of Hygiene & Tropical Medicine (AMR Centre)
   - Johns Hopkins Center for AMR Research
   - University of Oxford Big Data Institute

3. Research network participation:
   - DRIVE-AB consortium
   - GARDP surveillance network
   - Wellcome Trust AMR surveillance

PHASE 3: COMMERCIAL DATA (3-12 months)
-------------------------------------

1. IQVIA MIDAS Database:
   Contact: solutions@iqvia.com
   Budget: $75,000-200,000 annually
   Coverage: Global pharmaceutical sales

2. Alternative commercial sources:
   - Evaluate Pharma: Market intelligence
   - GlobalData Healthcare: Pharmaceutical market data
   - Pharma Intelligence: Regional market analysis

TECHNICAL IMPLEMENTATION
========================

1. Install required packages:
   pip install pandas openpyxl tabula-py requests beautifulsoup4

2. Create data directory structure:
   mkdir -p data/{ecdc,who,australia,national,commercial}

3. Set environment variables:
   export WHO_API_KEY="your_api_key"
   export IQVIA_USERNAME="your_username"  # If applicable
   export IQVIA_PASSWORD="your_password"  # If applicable

4. Run enhanced data collection:
   python enhanced_empirical_data_collection.py

EXPECTED OUTCOMES
================

After Phase 1: ~64 drug-region combinations with real data (22% coverage)
After Phase 2: ~150 drug-region combinations with real data (51% coverage)  
After Phase 3: ~260 drug-region combinations with real data (88% coverage)

Current coverage: ~3 combinations (1% coverage)
Target coverage: 260+ combinations (88%+ coverage)
"""