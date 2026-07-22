#!/usr/bin/env python3
"""
Observed-data ingestion and best-guess overlay generation for the AMR model.

Only records actually returned by an implemented external endpoint may be
classified as observed. Source-shaped generated values remain available as
explicit best-guess placeholders for diagnostics.

Configured source families:
1. WHO GLASS data (Global Antimicrobial Resistance and Use Surveillance System)
2. ECDC EARS-Net data (European Antimicrobial Resistance Surveillance Network)  
3. Australian NNDSS data (National Notifiable Diseases Surveillance System)
4. CDDEP ResistanceMap data (Center for Disease Dynamics, Economics & Policy)

Usage:
    python acquire_empirical_data.py --sources all
    python acquire_empirical_data.py --sources who,ecdc
    python acquire_empirical_data.py --sources who --years 2019-2023
"""

import requests
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import logging
from typing import List, Dict, Optional, Tuple
import json
from datetime import datetime
import time
import urllib.parse
import zipfile
import io
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EmpiricalDataAcquirer:
    """Ingest observations where implemented and label all generated fallbacks."""
    
    def __init__(self, base_dir: str = "data"):
        self.base_dir = Path(base_dir)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'AMR-Research-Tool/1.0 (Educational Research Use)'
        })
        
        # Create directory structure
        self.who_dir = self.base_dir / "who"
        self.ecdc_dir = self.base_dir / "ecdc" 
        self.australia_dir = self.base_dir / "australia"
        self.cddep_dir = self.base_dir / "cddep"
        
        for directory in [self.who_dir, self.ecdc_dir, self.australia_dir, self.cddep_dir]:
            directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _mark_generated_placeholder(
        df: pd.DataFrame,
        *,
        source_id: str,
        generation_method: str,
        rationale: str,
    ) -> pd.DataFrame:
        """Attach explicit non-observational provenance to generated values."""
        result = df.copy()
        result['overlay_provenance_class'] = 'generated_best_guess_placeholder'
        result['generated'] = True
        result['generation_method'] = generation_method
        result['source_id'] = source_id
        result['source_url_or_doi'] = ''
        result['reference_year'] = ''
        result['uncertainty'] = 'not_empirical_uncertainty'
        result['rationale'] = rationale
        result['last_reviewed'] = datetime.now().date().isoformat()
        return result

    @staticmethod
    def _mark_observed_download(
        df: pd.DataFrame,
        *,
        source_id: str,
        source_url: str,
        rationale: str,
    ) -> pd.DataFrame:
        """Attach the minimum provenance required for downloaded observations."""
        result = df.copy()
        result['overlay_provenance_class'] = 'observed_comparison'
        result['generated'] = False
        result['generation_method'] = 'not_generated_downloaded_api_data'
        result['source_id'] = source_id
        result['source_url_or_doi'] = source_url
        result['reference_year'] = (
            result['year'].astype(str) if 'year' in result.columns else ''
        )
        result['uncertainty'] = 'as_reported_or_not_available_from_endpoint'
        result['rationale'] = rationale
        result['last_reviewed'] = datetime.now().date().isoformat()
        return result
    
    def acquire_who_glass_data(self, years: List[int] = None) -> bool:
        """
        Acquire WHO GLASS surveillance data.
        
        WHO GLASS provides country-level AMR and AMU surveillance data.
        Access through WHO Global Health Observatory and GLASS dashboard.
        """
        logger.info("🏥 Acquiring WHO GLASS surveillance data...")
        
        if years is None:
            years = list(range(2017, 2024))  # GLASS started comprehensive reporting in 2017
            
        try:
            # Method 1: WHO Global Health Observatory (GHO) API
            success = self._acquire_who_gho_data(years)
            if success:
                logger.info("✅ WHO GLASS data acquisition completed")
                return True
                
            # Method 2: create source-shaped best-guess values (not observations)
            success = self._create_patterned_who_glass_data(years)
            if success:
                logger.info("WHO GLASS best-guess pattern dataset created")
                return True
                
            # Method 3: Create synthetic WHO GLASS structured data
            logger.warning("⚠️  Direct WHO GLASS access unavailable, creating structured data template")
            self._create_who_glass_template()
            return True
            
        except Exception as e:
            logger.error(f"❌ WHO GLASS acquisition failed: {e}")
            return False
    
    def _acquire_who_gho_data(self, years: List[int]) -> bool:
        """Acquire data from WHO Global Health Observatory."""
        
        # WHO GHO GLASS indicators
        indicators = [
            "GLASS_AMR_ECOLI_3GC",  # E. coli resistance to 3rd generation cephalosporins
            "GLASS_AMR_MRSA",       # MRSA resistance
            "GLASS_AMU_DDDPER1000", # Antimicrobial consumption DDD per 1000
        ]
        
        all_data = []
        
        for indicator in indicators:
            try:
                # WHO GHO REST API (documented public API)
                url = f"https://ghoapi.azureedge.net/api/{indicator}"
                
                logger.info(f"   Requesting {indicator}...")
                response = self.session.get(url, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    
                    # Process WHO GHO format
                    if 'value' in data:
                        for record in data['value']:
                            if record.get('TimeDim') and int(record['TimeDim']) in years:
                                processed_record = {
                                    'year': int(record['TimeDim']),
                                    'country': record.get('SpatialDim', ''),
                                    'indicator': indicator,
                                    'value': record.get('NumericValue'),
                                    'unit': record.get('UnitOfMeasure', ''),
                                    'source': 'WHO_GLASS_GHO'
                                }
                                all_data.append(processed_record)
                
                time.sleep(1)  # Rate limiting
                
            except Exception as e:
                logger.warning(f"   Failed to acquire {indicator}: {e}")
                continue
        
        if all_data:
            # Save WHO GLASS data
            df = self._mark_observed_download(
                pd.DataFrame(all_data),
                source_id='who_glass_gho_api',
                source_url='https://ghoapi.azureedge.net/api/',
                rationale='Records returned directly by the documented WHO GHO endpoint.',
            )
            output_path = self.who_dir / "glass_surveillance_data.csv"
            df.to_csv(output_path, index=False)
            logger.info(f"   Saved {len(df)} WHO GLASS records to {output_path}")
            return True
            
        return False
    
    def _create_patterned_who_glass_data(self, years: List[int]) -> bool:
        """Create source-shaped best-guess values; this does not access GLASS."""
        
        # This would access the Shiny app API endpoints if available
        # For now, create structured sample data based on GLASS reports
        
        logger.info("   Creating explicitly labelled WHO GLASS pattern placeholders...")
        
        # Sample GLASS-style data structure based on published reports
        glass_countries = [
            "Australia", "Austria", "Belgium", "Bulgaria", "Canada", "Croatia", 
            "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany",
            "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Lithuania",
            "Luxembourg", "Malta", "Netherlands", "Norway", "Poland", "Portugal",
            "Romania", "Slovak Republic", "Slovenia", "Spain", "Sweden", "United Kingdom"
        ]
        
        pathogens = ["Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus", 
                    "Streptococcus pneumoniae", "Enterococcus faecium", "Acinetobacter spp."]
        
        antibiotics = ["Ciprofloxacin", "Ceftriaxone", "Methicillin", "Vancomycin", 
                      "Carbapenem", "Colistin"]
        
        glass_data = []
        
        for year in years:
            for country in glass_countries:
                for pathogen in pathogens:
                    for antibiotic in antibiotics:
                        # Generate realistic resistance rates based on known patterns
                        if "Staphylococcus aureus" in pathogen and "Methicillin" in antibiotic:
                            # MRSA rates vary 1-50% across countries
                            resistance_rate = np.random.uniform(0.01, 0.50)
                        elif "Escherichia coli" in pathogen and "Ciprofloxacin" in antibiotic:
                            # E. coli fluoroquinolone resistance 5-35%
                            resistance_rate = np.random.uniform(0.05, 0.35)
                        elif "coli" in pathogen and ("Ceftriaxone" in antibiotic or "generation" in antibiotic):
                            # E. coli 3GC resistance 3-25%
                            resistance_rate = np.random.uniform(0.03, 0.25)
                        else:
                            resistance_rate = np.random.uniform(0.01, 0.30)
                        
                        glass_data.append({
                            'year': year,
                            'country': country,
                            'pathogen': pathogen,
                            'antibiotic': antibiotic,
                            'resistance_percentage': round(resistance_rate * 100, 2),
                            'number_tested': np.random.randint(50, 2000),
                            'source': 'WHO_GLASS',
                            'data_quality': 'representative' if np.random.random() > 0.3 else 'limited'
                        })
        
        # Save structured GLASS data
        df = self._mark_generated_placeholder(
            pd.DataFrame(glass_data),
            source_id='who_glass_pattern_template',
            generation_method='random_source_pattern_generator',
            rationale=(
                'Programmatically generated WHO-GLASS-structured values; the source '
                'name describes the intended format, not observational provenance.'
            ),
        )
        output_path = self.who_dir / "glass_amr_surveillance.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"   Created {len(df)} WHO GLASS-structured records at {output_path}")
        
        return True
    
    def _create_who_glass_template(self):
        """Create WHO GLASS data template with proper structure."""
        
        # Basic GLASS template based on official methodology
        template_data = {
            'country': ['Germany', 'France', 'Italy', 'Spain', 'Netherlands'],
            'year': [2022] * 5,
            'pathogen': ['Escherichia coli'] * 5,
            'antibiotic': ['Ciprofloxacin'] * 5,
            'resistance_percentage': [25.3, 28.7, 22.1, 31.2, 19.8],
            'number_tested': [1205, 987, 1456, 823, 1123],
            'source': ['WHO_GLASS'] * 5,
            'notes': ['Template data - replace with actual GLASS downloads'] * 5
        }
        
        df = self._mark_generated_placeholder(
            pd.DataFrame(template_data),
            source_id='who_glass_small_template',
            generation_method='fixed_source_pattern_template',
            rationale='Fixed illustrative template values; not downloaded observations.',
        )
        template_path = self.who_dir / "glass_template.csv"
        df.to_csv(template_path, index=False)
        logger.info(f"   Created WHO GLASS template at {template_path}")
    
    def acquire_ecdc_ears_net_data(self, years: List[int] = None) -> bool:
        """
        Acquire ECDC EARS-Net surveillance data.
        
        ECDC EARS-Net is the largest publicly-funded AMR surveillance in Europe.
        Data accessible through ECDC Surveillance Atlas.
        """
        logger.info("🇪🇺 Acquiring ECDC EARS-Net surveillance data...")
        
        if years is None:
            years = list(range(2018, 2024))
            
        try:
            # Method 1: ECDC Atlas API (not currently implemented)
            success = self._acquire_ecdc_atlas_data(years)
            if success:
                logger.info("✅ ECDC EARS-Net data acquisition completed")
                return True
                
            # Method 2: Annual report data extraction  
            success = self._acquire_ecdc_annual_reports(years)
            if success:
                logger.info("✅ ECDC annual report data acquisition completed") 
                return True
                
            # Method 3: create explicitly labelled source-patterned placeholders
            logger.warning("Direct ECDC observations unavailable; creating best-guess placeholders")
            self._create_comprehensive_ecdc_data(years)
            return True
            
        except Exception as e:
            logger.error(f"❌ ECDC acquisition failed: {e}")
            return False
    
    def _acquire_ecdc_atlas_data(self, years: List[int]) -> bool:
        """Access ECDC Surveillance Atlas API."""

        logger.info("   No ECDC Atlas API ingestion is implemented")
        return False
    
    def _create_comprehensive_ecdc_data(self, years: List[int]) -> bool:
        """Create comprehensive ECDC EARS-Net data based on published patterns."""
        
        # EU/EEA countries in EARS-Net
        eu_countries = [
            "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic",
            "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary",
            "Iceland", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta",
            "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Slovakia", 
            "Slovenia", "Spain", "Sweden"
        ]
        
        # EARS-Net priority pathogens
        ears_pathogens = [
            "Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
            "Acinetobacter species", "Streptococcus pneumoniae", "Enterococcus faecalis",
            "Enterococcus faecium", "Staphylococcus aureus"
        ]
        
        # Key antimicrobials per pathogen (EARS-Net standard)
        pathogen_antibiotics = {
            "Escherichia coli": ["Ampicillin", "Gentamicin", "Ciprofloxacin", "Ceftazidime", "Carbapenem"],
            "Klebsiella pneumoniae": ["Gentamicin", "Ciprofloxacin", "Ceftazidime", "Carbapenem", "Colistin"],
            "Pseudomonas aeruginosa": ["Gentamicin", "Ciprofloxacin", "Ceftazidime", "Carbapenem", "Colistin"],
            "Acinetobacter species": ["Gentamicin", "Ciprofloxacin", "Carbapenem", "Colistin"],
            "Streptococcus pneumoniae": ["Penicillin", "Erythromycin", "Levofloxacin"],
            "Enterococcus faecalis": ["Ampicillin", "Gentamicin", "Vancomycin"],
            "Enterococcus faecium": ["Ampicillin", "Vancomycin"],
            "Staphylococcus aureus": ["Methicillin", "Vancomycin"]
        }
        
        ecdc_data = []
        
        for year in years:
            for country in eu_countries:
                for pathogen in ears_pathogens:
                    for antibiotic in pathogen_antibiotics[pathogen]:
                        
                        # Generate realistic resistance rates based on ECDC trends
                        resistance_rate = self._get_realistic_resistance_rate(country, pathogen, antibiotic, year)
                        tested_isolates = np.random.randint(100, 5000)
                        
                        ecdc_data.append({
                            'year': year,
                            'country': country,
                            'pathogen': pathogen, 
                            'antibiotic': antibiotic,
                            'resistance_percentage': round(resistance_rate, 2),
                            'number_tested': tested_isolates,
                            'resistance_count': int(tested_isolates * resistance_rate / 100),
                            'source': 'ECDC_EARS_NET',
                            'specimen_type': 'Blood' if np.random.random() > 0.2 else 'CSF',
                            'reporting_level': 'National'
                        })
        
        # Save comprehensive ECDC data
        df = self._mark_generated_placeholder(
            pd.DataFrame(ecdc_data),
            source_id='ecdc_ears_net_pattern_template',
            generation_method='random_source_pattern_generator',
            rationale=(
                'Programmatically generated EARS-Net-structured values; the source '
                'name describes the intended format, not observational provenance.'
            ),
        )
        output_path = self.ecdc_dir / "ears_net_surveillance.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"   Created {len(df)} ECDC EARS-Net records at {output_path}")
        
        return True
    
    def _get_realistic_resistance_rate(self, country: str, pathogen: str, antibiotic: str, year: int) -> float:
        """Generate realistic resistance rates based on known epidemiological patterns."""
        
        # Base rates by pathogen-antibiotic combination (approximate European averages)
        base_rates = {
            ("Escherichia coli", "Ciprofloxacin"): 22.0,
            ("Escherichia coli", "Ceftazidime"): 12.0,
            ("Klebsiella pneumoniae", "Carbapenem"): 8.0,
            ("Staphylococcus aureus", "Methicillin"): 18.0,
            ("Enterococcus faecium", "Vancomycin"): 12.0,
            ("Pseudomonas aeruginosa", "Carbapenem"): 18.0,
            ("Acinetobacter species", "Carbapenem"): 35.0
        }
        
        # Get base rate
        key = (pathogen, antibiotic)
        base_rate = base_rates.get(key, 15.0)  # Default 15% if not specified
        
        # Country adjustments (simplified)
        country_factors = {
            "Germany": 0.8, "Netherlands": 0.6, "Denmark": 0.5, "Sweden": 0.5,
            "France": 1.0, "Italy": 1.3, "Spain": 1.2, "Poland": 1.4,
            "Romania": 1.8, "Bulgaria": 2.0, "Greece": 1.7
        }
        
        country_factor = country_factors.get(country, 1.0)
        
        # Year trend (slight increase over time)
        year_factor = 1.0 + (year - 2018) * 0.02
        
        # Add some random variation
        variation = np.random.uniform(0.8, 1.2)
        
        final_rate = base_rate * country_factor * year_factor * variation
        
        # Keep within realistic bounds
        return max(0.1, min(80.0, final_rate))
    
    def _acquire_ecdc_annual_reports(self, years: List[int]) -> bool:
        """Extract data from ECDC annual surveillance reports."""
        
        report_urls = {
            2023: "https://www.ecdc.europa.eu/en/publications-data/antimicrobial-resistance-eueea-ears-net-annual-epidemiological-report-2023",
            2022: "https://www.ecdc.europa.eu/en/publications-data/surveillance-antimicrobial-resistance-europe-2022",
            2021: "https://www.ecdc.europa.eu/en/publications-data/surveillance-antimicrobial-resistance-europe-2021"
        }
        
        for year in years:
            if year in report_urls:
                try:
                    logger.info(f"   Accessing ECDC {year} annual report...")
                    
                    # This would parse the actual PDF or web data
                    # For now, note the available report URLs
                    
                except Exception as e:
                    logger.warning(f"   Failed to process ECDC {year} report: {e}")
        
        return False  # Return False to trigger template creation
    
    def _create_ecdc_template(self):
        """Create ECDC template with current sample data expanded."""
        
        # Expand the existing minimal ECDC data
        current_data = {
            'Country': ['Germany', 'France', 'Spain', 'Germany', 'Netherlands'],
            'Year': [2022, 2022, 2022, 2021, 2022],  
            'Bacteria': ['Escherichia coli', 'Escherichia coli', 'E. coli', 'S. aureus', 'K. pneumoniae'],
            'Antibiotic': ['Ciprofloxacin', 'Ciprofloxacin', 'Tetracycline', 'Methicillin', 'Carbapenem'],
            'Resistance_percentage': [25.3, 28.7, 35.2, 8.1, 12.4],
            'Number_tested': [1205, 987, 756, 2134, 1567]
        }
        
        df = self._mark_generated_placeholder(
            pd.DataFrame(current_data),
            source_id='ecdc_ears_net_small_template',
            generation_method='fixed_source_pattern_template',
            rationale='Fixed illustrative template values; not downloaded observations.',
        )
        template_path = self.ecdc_dir / "ecdc_template_expanded.csv"
        df.to_csv(template_path, index=False)
        logger.info(f"   Created expanded ECDC template at {template_path}")
    
    def acquire_australian_nndss_data(self, years: List[int] = None) -> bool:
        """
        Acquire Australian NNDSS surveillance data.
        
        National Notifiable Diseases Surveillance System data for AMR.
        """
        logger.info("🇦🇺 Acquiring Australian NNDSS surveillance data...")
        
        if years is None:
            years = list(range(2019, 2024))
            
        try:
            # No observational NNDSS ingestion is implemented.
            success = self._create_patterned_nndss_data(years)
            if success:
                logger.info("Australian NNDSS best-guess pattern dataset created")
                return True
                
            # Create Australian surveillance template
            logger.warning("⚠️  Direct NNDSS access unavailable, creating template")
            self._create_australian_template()
            return True
            
        except Exception as e:
            logger.error(f"❌ Australian NNDSS acquisition failed: {e}")
            return False
    
    def _create_patterned_nndss_data(self, years: List[int]) -> bool:
        """Create source-shaped best-guess values; this does not access NNDSS."""
        
        # Australian surveillance data structure
        australian_states = ["NSW", "VIC", "QLD", "WA", "SA", "TAS", "ACT", "NT"]
        
        # Key pathogens under Australian surveillance
        aus_pathogens = [
            "Staphylococcus aureus", "Enterococcus species", "Escherichia coli",
            "Klebsiella pneumoniae", "Pseudomonas aeruginosa"
        ]
        
        aus_data = []
        
        for year in years:
            for state in australian_states:
                for pathogen in aus_pathogens:
                    # Generate Australian-specific resistance patterns
                    if "aureus" in pathogen:  # MRSA rates in Australia
                        resistance_rate = np.random.uniform(15.0, 25.0)
                        antibiotic = "Methicillin"
                    elif "coli" in pathogen:  # E. coli fluoroquinolone resistance
                        resistance_rate = np.random.uniform(8.0, 18.0)
                        antibiotic = "Ciprofloxacin"
                    else:
                        resistance_rate = np.random.uniform(5.0, 20.0)
                        antibiotic = "Various"
                    
                    aus_data.append({
                        'year': year,
                        'state': state,
                        'pathogen': pathogen,
                        'antibiotic': antibiotic,
                        'resistance_percentage': round(resistance_rate, 2),
                        'number_tested': np.random.randint(50, 800),
                        'source': 'NNDSS_Australia',
                        'reporting_level': 'State'
                    })
        
        # Save Australian data
        df = self._mark_generated_placeholder(
            pd.DataFrame(aus_data),
            source_id='australian_nndss_pattern_template',
            generation_method='random_source_pattern_generator',
            rationale=(
                'Programmatically generated NNDSS-structured values; the source '
                'name describes the intended format, not observational provenance.'
            ),
        )
        output_path = self.australia_dir / "nndss_surveillance.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"   Created {len(df)} Australian NNDSS records at {output_path}")
        
        return True
    
    def _create_australian_template(self):
        """Create Australian surveillance data template."""
        
        template_data = {
            'year': [2022] * 5,
            'state': ['NSW', 'VIC', 'QLD', 'WA', 'SA'],
            'pathogen': ['Staphylococcus aureus'] * 5,
            'antibiotic': ['Methicillin'] * 5,
            'resistance_percentage': [18.5, 20.1, 17.8, 19.6, 21.2],
            'number_tested': [456, 623, 389, 512, 334],
            'source': ['NNDSS'] * 5,
            'notes': ['Template - replace with NNDSS downloads'] * 5
        }
        
        df = self._mark_generated_placeholder(
            pd.DataFrame(template_data),
            source_id='australian_nndss_small_template',
            generation_method='fixed_source_pattern_template',
            rationale='Fixed illustrative template values; not downloaded observations.',
        )
        template_path = self.australia_dir / "nndss_template.csv" 
        df.to_csv(template_path, index=False)
        logger.info(f"   Created Australian template at {template_path}")
    
    def acquire_cddep_resistancemap_data(self, years: List[int] = None) -> bool:
        """
        Acquire CDDEP ResistanceMap data.
        
        Center for Disease Dynamics, Economics & Policy global resistance data.
        """
        logger.info("🌍 Acquiring CDDEP ResistanceMap data...")
        
        if years is None:
            years = list(range(2018, 2024))
            
        try:
            # CDDEP ResistanceMap API/portal
            success = self._create_patterned_resistancemap_data(years)
            if success:
                logger.info("ResistanceMap best-guess pattern dataset created")
                return True
                
            # Create global resistance mapping template
            logger.warning("⚠️  Direct ResistanceMap access unavailable, creating global template")
            self._create_cddep_template()
            return True
            
        except Exception as e:
            logger.error(f"❌ CDDEP acquisition failed: {e}")
            return False
    
    def _create_patterned_resistancemap_data(self, years: List[int]) -> bool:
        """Create source-shaped best-guess values; this does not access ResistanceMap."""
        
        # Global countries with AMR data
        global_countries = [
            "United States", "Canada", "Brazil", "Argentina", "United Kingdom", "Germany", 
            "France", "Italy", "Spain", "Netherlands", "Russia", "Turkey", "India", 
            "China", "Japan", "South Korea", "Australia", "South Africa", "Kenya", "Egypt"
        ]
        
        # Global resistance patterns
        global_data = []
        
        for year in years:
            for country in global_countries:
                # Key global pathogens
                for pathogen in ["E. coli", "K. pneumoniae", "S. aureus", "S. pneumoniae"]:
                    for antibiotic in ["Fluoroquinolones", "3rd-gen Cephalosporins", "Methicillin", "Penicillin"]:
                        
                        # Global resistance patterns vary significantly by region
                        base_rate = self._get_global_resistance_rate(country, pathogen, antibiotic)
                        
                        global_data.append({
                            'year': year,
                            'country': country,
                            'pathogen': pathogen,
                            'antibiotic': antibiotic,
                            'resistance_percentage': round(base_rate, 2),
                            'number_tested': np.random.randint(100, 2000),
                            'source': 'CDDEP_ResistanceMap',
                            'region': self._get_region(country),
                            'income_level': self._get_income_level(country)
                        })
        
        # Save CDDEP data
        df = self._mark_generated_placeholder(
            pd.DataFrame(global_data),
            source_id='cddep_resistancemap_pattern_template',
            generation_method='random_source_pattern_generator',
            rationale=(
                'Programmatically generated ResistanceMap-structured values; the '
                'source name describes the intended format, not observational provenance.'
            ),
        )
        output_path = self.cddep_dir / "resistancemap_surveillance.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"   Created {len(df)} CDDEP ResistanceMap records at {output_path}")
        
        return True
    
    def _get_global_resistance_rate(self, country: str, pathogen: str, antibiotic: str) -> float:
        """Generate realistic global resistance rates by region."""
        
        # Regional base rates (simplified)
        regional_rates = {
            "High-income": {"E. coli": 15.0, "K. pneumoniae": 12.0, "S. aureus": 20.0},
            "Upper-middle": {"E. coli": 25.0, "K. pneumoniae": 20.0, "S. aureus": 30.0}, 
            "Lower-middle": {"E. coli": 35.0, "K. pneumoniae": 30.0, "S. aureus": 40.0},
            "Low-income": {"E. coli": 45.0, "K. pneumoniae": 40.0, "S. aureus": 50.0}
        }
        
        income_level = self._get_income_level(country)
        base_rate = regional_rates[income_level].get(pathogen, 25.0)
        
        # Antibiotic-specific adjustments
        if "Methicillin" in antibiotic and "aureus" in pathogen:
            base_rate *= 0.9  # MRSA rates
        elif "Fluoroquinolone" in antibiotic:
            base_rate *= 1.2  # Generally higher FQ resistance
        elif "Cephalosporin" in antibiotic:
            base_rate *= 0.8  # Moderate cephalosporin resistance
            
        # Add variation
        return max(1.0, min(80.0, base_rate * np.random.uniform(0.7, 1.3)))
    
    def _get_region(self, country: str) -> str:
        """Get WHO region for country."""
        regions = {
            "Americas": ["United States", "Canada", "Brazil", "Argentina"],
            "Europe": ["United Kingdom", "Germany", "France", "Italy", "Spain", "Netherlands", "Russia"],
            "South-East Asia": ["India"],
            "Western Pacific": ["China", "Japan", "South Korea", "Australia"], 
            "Africa": ["South Africa", "Kenya"],
            "Eastern Mediterranean": ["Turkey", "Egypt"]
        }
        
        for region, countries in regions.items():
            if country in countries:
                return region
        return "Other"
    
    def _get_income_level(self, country: str) -> str:
        """Get World Bank income classification."""
        income_levels = {
            "High-income": ["United States", "Canada", "United Kingdom", "Germany", "France", 
                           "Italy", "Spain", "Netherlands", "Japan", "South Korea", "Australia"],
            "Upper-middle": ["Brazil", "Argentina", "Russia", "Turkey", "China", "South Africa"],
            "Lower-middle": ["India", "Egypt", "Kenya"],
            "Low-income": []  # Add as needed
        }
        
        for level, countries in income_levels.items():
            if country in countries:
                return level
        return "Upper-middle"  # Default
    
    def _create_cddep_template(self):
        """Create CDDEP global resistance mapping template."""
        
        template_data = {
            'year': [2022] * 8,
            'country': ['India', 'China', 'United States', 'Brazil', 'Germany', 'South Africa', 'Kenya', 'Turkey'],
            'pathogen': ['E. coli'] * 8,
            'antibiotic': ['Fluoroquinolones'] * 8,
            'resistance_percentage': [65.2, 58.7, 22.1, 45.8, 18.9, 72.3, 78.5, 52.1],
            'number_tested': [856, 1245, 2156, 623, 1567, 345, 234, 789],
            'source': ['CDDEP_ResistanceMap'] * 8,
            'region': ['South-East Asia', 'Western Pacific', 'Americas', 'Americas', 'Europe', 'Africa', 'Africa', 'Eastern Mediterranean'],
            'notes': ['Template - replace with ResistanceMap downloads'] * 8
        }
        
        df = self._mark_generated_placeholder(
            pd.DataFrame(template_data),
            source_id='cddep_resistancemap_small_template',
            generation_method='fixed_source_pattern_template',
            rationale='Fixed illustrative template values; not downloaded observations.',
        )
        template_path = self.cddep_dir / "resistancemap_template.csv"
        df.to_csv(template_path, index=False) 
        logger.info(f"   Created CDDEP template at {template_path}")
    
    def generate_integration_summary(self) -> None:
        """Generate summary report of acquired empirical data."""
        
        logger.info("📊 Generating empirical data integration summary...")
        
        summary = {
            'summary_contract_version': 1,
            'summary_date': datetime.now().isoformat(),
            'sources_present': [],
            'total_records': 0,
            'observed_records': 0,
            'best_guess_placeholder_records': 0,
            'unclassified_records': 0,
            'coverage_assessment': {}
        }
        
        # Check each data source
        data_sources = [
            ('WHO GLASS', self.who_dir),
            ('ECDC EARS-Net', self.ecdc_dir), 
            ('Australian NNDSS', self.australia_dir),
            ('CDDEP ResistanceMap', self.cddep_dir)
        ]
        known_generated_files = {
            'glass_amr_surveillance.csv',
            'glass_template.csv',
            'ears_net_surveillance.csv',
            'ecdc_template_expanded.csv',
            'nndss_surveillance.csv',
            'nndss_template.csv',
            'resistancemap_surveillance.csv',
            'resistancemap_template.csv',
        }
        
        for source_name, source_dir in data_sources:
            csv_files = list(source_dir.glob("*.csv"))
            
            if csv_files:
                total_source_records = 0
                observed_records = 0
                placeholder_records = 0
                unclassified_records = 0
                for csv_file in csv_files:
                    try:
                        df = pd.read_csv(csv_file)
                        total_source_records += len(df)
                        if 'overlay_provenance_class' in df.columns:
                            classes = df['overlay_provenance_class'].fillna('')
                            observed_records += int((classes == 'observed_comparison').sum())
                            placeholder_records += int(classes.isin([
                                'generated_best_guess_placeholder',
                                'source_informed_best_guess_placeholder_provenance_unverified'
                            ]).sum())
                            unclassified_records += int((classes == '').sum())
                        else:
                            if csv_file.name in known_generated_files:
                                placeholder_records += len(df)
                            else:
                                unclassified_records += len(df)
                    except Exception as e:
                        logger.warning(f"Could not read {csv_file}: {e}")
                
                summary['sources_present'].append(source_name)
                summary['total_records'] += total_source_records
                summary['observed_records'] += observed_records
                summary['best_guess_placeholder_records'] += placeholder_records
                summary['unclassified_records'] += unclassified_records
                if observed_records and placeholder_records:
                    status = 'mixed_observed_and_placeholder'
                elif observed_records:
                    status = 'observed'
                elif placeholder_records:
                    status = 'best_guess_placeholders_only'
                else:
                    status = 'provenance_unclassified'
                summary['coverage_assessment'][source_name] = {
                    'files': len(csv_files),
                    'records': total_source_records,
                    'observed_records': observed_records,
                    'best_guess_placeholder_records': placeholder_records,
                    'unclassified_records': unclassified_records,
                    'status': status
                }
                
                logger.info(f"   ✅ {source_name}: {total_source_records:,} records in {len(csv_files)} files")
            else:
                summary['coverage_assessment'][source_name] = {
                    'files': 0, 
                    'records': 0,
                    'status': 'not_acquired'
                }
                logger.info(f"   ❌ {source_name}: No data acquired")
        
        # Save summary
        summary_path = self.base_dir / "empirical_data_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"📋 Empirical data summary saved to {summary_path}")
        logger.info(
            "Overlay records: %s observed, %s best-guess placeholders, %s unclassified",
            f"{summary['observed_records']:,}",
            f"{summary['best_guess_placeholder_records']:,}",
            f"{summary['unclassified_records']:,}",
        )
        
        return summary

def main():
    """Main execution function with command-line interface."""
    
    parser = argparse.ArgumentParser(description="Acquire empirical AMR surveillance data")
    parser.add_argument("--sources", default="all", 
                       help="Data sources to acquire: all,who,ecdc,australia,cddep")
    parser.add_argument("--years", default="2019-2023",
                       help="Year range: 2019-2023 or specific years: 2021,2022,2023")
    parser.add_argument("--output-dir", default="data",
                       help="Output directory for acquired data")
    
    args = parser.parse_args()
    
    # Parse years
    if '-' in args.years:
        start_year, end_year = map(int, args.years.split('-'))
        years = list(range(start_year, end_year + 1))
    else:
        years = [int(y.strip()) for y in args.years.split(',')]
    
    # Parse sources 
    if args.sources.lower() == 'all':
        sources = ['who', 'ecdc', 'australia', 'cddep']
    else:
        sources = [s.strip().lower() for s in args.sources.split(',')]
    
    # Initialize acquirer
    acquirer = EmpiricalDataAcquirer(args.output_dir)
    
    logger.info(f"🚀 Starting empirical data acquisition for {sources} ({years})")
    
    # Acquire each requested source
    success_count = 0
    
    if 'who' in sources:
        if acquirer.acquire_who_glass_data(years):
            success_count += 1
    
    if 'ecdc' in sources:
        if acquirer.acquire_ecdc_ears_net_data(years):
            success_count += 1
            
    if 'australia' in sources:
        if acquirer.acquire_australian_nndss_data(years):
            success_count += 1
            
    if 'cddep' in sources:
        if acquirer.acquire_cddep_resistancemap_data(years):
            success_count += 1
    
    # Generate summary
    summary = acquirer.generate_integration_summary()
    
    logger.info(f"🏁 Data acquisition completed: {success_count}/{len(sources)} sources successful")
    
    return success_count == len(sources)

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
