#!/usr/bin/env python3
"""
Enhanced Empirical Data Sourcing Strategy
Comprehensive plan for expanding drug consumption data coverage across all regions and drugs
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
import requests
import json

class EnhancedDataStrategy:
    """
    Strategy for sourcing real empirical data for all 49 drugs across 6 regions
    """
    
    def __init__(self):
        # Your 49 drugs from the simulation
        self.drugs = [
            'sulfanilamide', 'penicilling', 'ampicillin', 'amoxicillin',
            'piperacillin', 'ticarcillin', 'cephalexin', 'cefazolin',
            'cefuroxime', 'ceftriaxone', 'ceftazidime', 'cefepime', 'ceftaroline',
            'meropenem', 'imipenem_c', 'ertapenem', 'aztreonam', 'erythromycin',
            'azithromycin', 'clarithromycin', 'clindamycin', 'gentamicin',
            'tobramycin', 'amikacin', 'ciprofloxacin', 'levofloxacin',
            'moxifloxacin', 'ofloxacin', 'tetracycline', 'doxyclycline',
            'minocycline', 'vancomycin', 'teicoplanin', 'linezolid',
            'tedizolid', 'quinu_dalfo', 'trim_sulf', 'chlorampheni',
            'nitrofurantoin', 'retapamulin', 'fusidic_a', 'metronidazole',
            'furazolidone', 'rifampicin', 'amoxicillin_clavulanate',
            'piperacillin_tazobactam', 'ampicillin_sulbactam',
            'ticarcillin_clavulanate', 'ceftazidime_avibactam',
            'meropenem_vaborbactam', 'colistin'
        ]
        
        self.regions = ['north_america', 'europe', 'asia', 'africa', 'south_america', 'oceania']
        
    def get_comprehensive_data_sources(self) -> Dict:
        """
        Complete mapping of data sources for each region and drug coverage
        """
        
        data_sources = {
            # === EUROPEAN REGION ===
            'europe': {
                'primary_sources': {
                    'ecdc_esac_net': {
                        'description': 'European Surveillance of Antimicrobial Consumption Network',
                        'url': 'https://www.ecdc.europa.eu/en/antimicrobial-consumption/database',
                        'coverage': '31 EU/EEA countries',
                        'data_type': 'DDD per 1000 inhabitants per day',
                        'drugs_covered': [
                            'penicilling', 'ampicillin', 'amoxicillin', 'amoxicillin_clavulanate',
                            'cephalexin', 'cefuroxime', 'ceftriaxone', 'erythromycin',
                            'azithromycin', 'clarithromycin', 'clindamycin', 'ciprofloxacin',
                            'levofloxacin', 'moxifloxacin', 'tetracycline', 'doxyclycline',
                            'trim_sulf', 'nitrofurantoin', 'vancomycin', 'linezolid',
                            'metronidazole', 'rifampicin'
                        ],
                        'years_available': '2012-2023',
                        'access': 'Free public download',
                        'file_format': 'Excel/CSV'
                    },
                    'european_medicines_agency': {
                        'description': 'EMA consumption reporting',
                        'url': 'https://www.ema.europa.eu/en/veterinary-regulatory/overview/antimicrobial-resistance',
                        'coverage': 'EU member states',
                        'drugs_covered': ['most major antibiotics'],
                        'access': 'Public reports'
                    }
                },
                'secondary_sources': {
                    'national_agencies': {
                        'germany_dimdi': 'https://www.dimdi.de/dynamic/de/arzneimittel/arzneimittel-informationssystem/',
                        'france_ansm': 'https://ansm.sante.fr/page/les-donnees-de-consommation-et-de-resistance-aux-antibiotiques-en-france',
                        'uk_ons': 'https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare',
                        'netherlands_swab': 'https://swab.nl/nl/',
                        'sweden_strama': 'https://strama.se/'
                    }
                }
            },
            
            # === NORTH AMERICA ===
            'north_america': {
                'primary_sources': {
                    'iqvia_midas': {
                        'description': 'IQVIA MIDAS pharmaceutical sales database',
                        'coverage': 'USA, Canada, Mexico',
                        'drugs_covered': 'All major antibiotics',
                        'access': 'Commercial license required',
                        'cost': '$50,000-200,000 annually',
                        'contact': 'https://www.iqvia.com/solutions/commercialization/brand-strategy-and-management/midas'
                    },
                    'cdc_narms': {
                        'description': 'CDC National Antimicrobial Resistance Monitoring System',
                        'url': 'https://www.cdc.gov/narms/index.html',
                        'coverage': 'United States',
                        'access': 'Free public data'
                    },
                    'canada_phac': {
                        'description': 'Public Health Agency of Canada',
                        'url': 'https://www.canada.ca/en/public-health/services/antibiotic-antimicrobial-resistance.html',
                        'coverage': 'Canada',
                        'access': 'Public reports'
                    }
                },
                'alternative_sources': {
                    'prescription_databases': {
                        'us_medicare': 'https://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data/Part-D-Prescriber.html',
                        'us_medicaid': 'https://www.medicaid.gov/medicaid/prescription-drugs/index.html'
                    }
                }
            },
            
            # === ASIA ===
            'asia': {
                'primary_sources': {
                    'japan_niid': {
                        'description': 'National Institute of Infectious Diseases Japan',
                        'url': 'https://www.niid.go.jp/niid/en/antimicrobial-resistance-en.html',
                        'coverage': 'Japan',
                        'drugs_covered': 'Major antibiotics',
                        'access': 'Public surveillance reports'
                    },
                    'south_korea_kdca': {
                        'description': 'Korea Disease Control and Prevention Agency',
                        'url': 'https://www.kdca.go.kr/contents.es?mid=a20301000000',
                        'coverage': 'South Korea',
                        'access': 'Government reports'
                    },
                    'china_nhc': {
                        'description': 'National Health Commission of China',
                        'coverage': 'China',
                        'note': 'Limited public access, may require academic collaboration'
                    },
                    'who_searo': {
                        'description': 'WHO South-East Asia Regional Office',
                        'url': 'https://www.who.int/southeastasia/health-topics/antimicrobial-resistance',
                        'coverage': 'Bangladesh, India, Indonesia, Thailand, etc.',
                        'access': 'Regional surveillance reports'
                    }
                },
                'proxy_sources': {
                    'pharmaceutical_industry': {
                        'indian_pharma_exports': 'https://www.ibef.org/industry/pharmaceutical-india.aspx',
                        'china_nmpa': 'https://www.nmpa.gov.cn/',
                        'asean_pharma': 'Various national regulatory authorities'
                    }
                }
            },
            
            # === AFRICA ===
            'africa': {
                'primary_sources': {
                    'who_afro': {
                        'description': 'WHO Regional Office for Africa',
                        'url': 'https://www.afro.who.int/health-topics/antimicrobial-resistance',
                        'coverage': '47 African countries',
                        'access': 'Regional surveillance reports'
                    },
                    'africa_cdc': {
                        'description': 'Africa Centres for Disease Control and Prevention',
                        'url': 'https://africacdc.org/programme/laboratory-systems-and-networks/',
                        'coverage': 'Continental surveillance'
                    },
                    'national_programs': {
                        'south_africa_nicd': 'https://www.nicd.ac.za/centres/centre-for-healthcare-associated-infections-antimicrobial-resistance-and-mycoses/',
                        'nigeria_ncdc': 'https://ncdc.gov.ng/',
                        'kenya_kemri': 'https://www.kemri.org/',
                        'ghana_ghs': 'https://www.ghanahealthservice.org/'
                    }
                },
                'proxy_sources': {
                    'pharmaceutical_imports': {
                        'description': 'Use import/export statistics as proxy for consumption',
                        'sources': ['National customs data', 'WHO pharmaceutical imports database']
                    }
                }
            },
            
            # === SOUTH AMERICA ===
            'south_america': {
                'primary_sources': {
                    'paho': {
                        'description': 'Pan American Health Organization',
                        'url': 'https://www.paho.org/en/topics/antimicrobial-resistance',
                        'coverage': 'Latin America and Caribbean',
                        'access': 'Regional surveillance reports'
                    },
                    'national_programs': {
                        'brazil_anvisa': 'https://www.gov.br/anvisa/pt-br/assuntos/medicamentos/controle-de-medicamentos',
                        'argentina_anmat': 'https://www.argentina.gob.ar/anmat',
                        'colombia_invima': 'https://www.invima.gov.co/',
                        'chile_isp': 'https://www.ispch.cl/'
                    }
                },
                'regional_networks': {
                    'relavra': {
                        'description': 'Latin American Network for Antimicrobial Resistance Surveillance',
                        'coverage': 'Multiple Latin American countries'
                    }
                }
            },
            
            # === OCEANIA ===
            'oceania': {
                'primary_sources': {
                    'australia_aura': {
                        'description': 'Antimicrobial Use and Resistance in Australia (AURA)',
                        'url': 'https://www.safetyandquality.gov.au/our-work/antimicrobial-resistance/antimicrobial-use-and-resistance-australia-surveillance-system',
                        'coverage': 'Australia',
                        'drugs_covered': 'Comprehensive antibiotic coverage',
                        'access': 'Public surveillance reports'
                    },
                    'new_zealand_esr': {
                        'description': 'Institute of Environmental Science and Research',
                        'url': 'https://www.esr.cri.nz/our-expertise/public-health/',
                        'coverage': 'New Zealand',
                        'access': 'Public health surveillance'
                    }
                }
            }
        }
        
        return data_sources
    
    def create_data_acquisition_plan(self) -> Dict:
        """
        Step-by-step plan for acquiring empirical data
        """
        
        plan = {
            'phase_1_immediate': {
                'description': 'Expand existing free/public sources',
                'timeline': '1-2 months',
                'actions': [
                    {
                        'task': 'Download complete ECDC ESAC-Net database',
                        'url': 'https://www.ecdc.europa.eu/en/antimicrobial-consumption/database',
                        'expected_drugs': 22,
                        'expected_countries': 31,
                        'years': '2012-2023'
                    },
                    {
                        'task': 'Collect WHO regional surveillance reports',
                        'sources': ['WHO AFRO', 'WHO SEARO', 'WHO PAHO'],
                        'expected_coverage': 'Regional aggregates for major antibiotics'
                    },
                    {
                        'task': 'Download Australia AURA surveillance data',
                        'url': 'https://www.safetyandquality.gov.au/our-work/antimicrobial-resistance/antimicrobial-use-and-resistance-australia-surveillance-system',
                        'expected_drugs': 15-20,
                        'coverage': 'Australia complete'
                    },
                    {
                        'task': 'Collect national surveillance reports',
                        'countries': ['Japan', 'South Korea', 'Canada', 'Brazil', 'South Africa'],
                        'format': 'Manual extraction from PDFs'
                    }
                ]
            },
            
            'phase_2_commercial': {
                'description': 'Acquire commercial pharmaceutical databases',
                'timeline': '3-6 months',
                'budget_required': '$50,000-100,000',
                'actions': [
                    {
                        'task': 'IQVIA MIDAS global expansion',
                        'coverage': 'All regions, all major antibiotics',
                        'cost': '$50,000-200,000 annually',
                        'contact': 'IQVIA sales team'
                    },
                    {
                        'task': 'GERS (Global Epidemiology Reporting Services)',
                        'description': 'Alternative to IQVIA',
                        'coverage': 'Global pharmaceutical sales'
                    },
                    {
                        'task': 'Pharmaprojects database',
                        'description': 'Drug development and sales tracking',
                        'focus': 'Pipeline antibiotics'
                    }
                ]
            },
            
            'phase_3_partnerships': {
                'description': 'Academic and institutional collaborations',
                'timeline': '6-12 months',
                'actions': [
                    {
                        'task': 'WHO collaboration agreement',
                        'benefit': 'Access to unpublished surveillance data',
                        'contact': 'WHO Global Action Plan on AMR focal points'
                    },
                    {
                        'task': 'Academic partnerships',
                        'institutions': [
                            'London School of Hygiene & Tropical Medicine',
                            'Johns Hopkins Bloomberg School of Public Health',
                            'University of Oxford - Big Data Institute',
                            'Wellcome Trust AMR surveillance network'
                        ],
                        'benefit': 'Research data sharing agreements'
                    },
                    {
                        'task': 'Pharmaceutical industry partnerships',
                        'companies': ['Pfizer', 'GSK', 'Roche', 'Novartis'],
                        'benefit': 'Sales data for specific compounds'
                    }
                ]
            }
        }
        
        return plan
    
    def estimate_data_coverage_improvement(self) -> Dict:
        """
        Estimate how much empirical coverage we can achieve
        """
        
        coverage_estimates = {
            'current_state': {
                'total_drug_region_combinations': len(self.drugs) * len(self.regions),  # 49 * 6 = 294
                'empirical_coverage': 3,  # Only ~3 combinations have real data
                'percentage_empirical': 1.0
            },
            
            'after_phase_1': {
                'europe_expansion': {
                    'ecdc_drugs': 22,  # Major antibiotics covered by ECDC
                    'countries': 31,
                    'regional_coverage': 22  # 22 drugs for Europe region
                },
                'oceania_expansion': {
                    'australia_drugs': 18,  # AURA surveillance
                    'regional_coverage': 18
                },
                'partial_other_regions': {
                    'who_reports': 8,  # Major antibiotics in regional reports
                    'regions_covered': 3,  # Asia, Africa, South America
                    'combinations': 24  # 8 drugs * 3 regions
                },
                'total_empirical': 22 + 18 + 24,  # 64 combinations
                'percentage_empirical': 21.8
            },
            
            'after_phase_2': {
                'iqvia_expansion': {
                    'drugs_covered': 35,  # Most antibiotics have sales data
                    'regions_covered': 6,  # Global coverage
                    'combinations': 210  # 35 * 6
                },
                'total_empirical': 210,
                'percentage_empirical': 71.4
            },
            
            'after_phase_3': {
                'academic_partnerships': {
                    'additional_combinations': 50,  # Fill remaining gaps
                    'total_empirical': 260,
                    'percentage_empirical': 88.4
                }
            }
        }
        
        return coverage_estimates

def main():
    """
    Generate comprehensive data sourcing strategy
    """
    
    strategy = EnhancedDataStrategy()
    
    print("🌍 ENHANCED EMPIRICAL DATA SOURCING STRATEGY")
    print("=" * 70)
    
    # Show data sources
    sources = strategy.get_comprehensive_data_sources()
    print(f"\n📊 DATA SOURCES BY REGION:")
    for region, data in sources.items():
        print(f"\n🌍 {region.upper().replace('_', ' ')}:")
        for source_type, source_info in data.items():
            print(f"  {source_type}:")
            if isinstance(source_info, dict):
                for source, details in source_info.items():
                    if isinstance(details, dict):
                        print(f"    • {source}: {details.get('description', 'Data source')}")
                    else:
                        print(f"    • {source}: {details}")
    
    # Show acquisition plan
    plan = strategy.create_data_acquisition_plan()
    print(f"\n📋 DATA ACQUISITION PLAN:")
    for phase, details in plan.items():
        print(f"\n🎯 {phase.upper().replace('_', ' ')}:")
        print(f"   Timeline: {details['timeline']}")
        if 'budget_required' in details:
            print(f"   Budget: {details['budget_required']}")
        print(f"   Actions:")
        for action in details['actions']:
            print(f"     • {action['task']}")
    
    # Show coverage improvements
    coverage = strategy.estimate_data_coverage_improvement()
    print(f"\n📈 EXPECTED COVERAGE IMPROVEMENTS:")
    for phase, data in coverage.items():
        if 'percentage_empirical' in data:
            print(f"  {phase}: {data['percentage_empirical']:.1f}% empirical coverage")

if __name__ == "__main__":
    main()