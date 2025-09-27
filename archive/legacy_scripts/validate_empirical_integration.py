#!/usr/bin/env python3
"""
Empirical Data Integration Validation Summary

This script validates the successful integration of enhanced empirical data sources
into the AMR simulation analysis system.
"""

import pandas as pd
from pathlib import Path
import json
from datetime import datetime

def validate_empirical_integration():
    """Validate the enhanced empirical data integration."""
    
    print("🔍 Validating Enhanced Empirical Data Integration")
    print("=" * 60)
    
    validation_results = {
        'validation_timestamp': datetime.now().isoformat(),
        'data_sources_validated': [],
        'integration_status': {},
        'plot_generation_status': {},
        'data_quality_metrics': {},
        'coverage_assessment': {}
    }
    
    # 1. Validate Data Acquisition
    print("\n1. Data Acquisition Validation")
    print("-" * 40)
    
    data_sources = [
        ('WHO GLASS', Path('data/who/glass_amr_surveillance.csv')),
        ('ECDC EARS-Net', Path('data/ecdc/ears_net_surveillance.csv')),
        ('Australian NNDSS', Path('data/australia/nndss_surveillance.csv')),
        ('CDDEP ResistanceMap', Path('data/cddep/resistancemap_surveillance.csv'))
    ]
    
    total_surveillance_records = 0
    
    for source_name, file_path in data_sources:
        if file_path.exists():
            try:
                df = pd.read_csv(file_path)
                record_count = len(df)
                total_surveillance_records += record_count
                
                validation_results['data_sources_validated'].append(source_name)
                validation_results['integration_status'][source_name] = {
                    'status': 'success',
                    'records': record_count,
                    'file_path': str(file_path),
                    'columns': list(df.columns)
                }
                
                print(f"✅ {source_name:20} | {record_count:6,} records | {file_path}")
                
                # Sample data quality check
                if 'year' in df.columns and 'resistance_percentage' in df.columns:
                    year_range = f"{df['year'].min()}-{df['year'].max()}"
                    resistance_range = f"{df['resistance_percentage'].min():.1f}-{df['resistance_percentage'].max():.1f}%"
                    print(f"   📊 Years: {year_range} | Resistance: {resistance_range}")
                    
            except Exception as e:
                validation_results['integration_status'][source_name] = {
                    'status': 'error',
                    'error': str(e)
                }
                print(f"❌ {source_name:20} | Error: {e}")
        else:
            validation_results['integration_status'][source_name] = {
                'status': 'missing',
                'file_path': str(file_path)
            }
            print(f"⚠️  {source_name:20} | File not found: {file_path}")
    
    print(f"\n📈 Total surveillance records acquired: {total_surveillance_records:,}")
    
    # 2. Validate Enhanced Loader Integration
    print("\n2. Enhanced Loader Integration Validation")
    print("-" * 40)
    
    try:
        # Test the enhanced loader
        from enhanced_empirical_loader import load_integrated_empirical_data
        
        integrated_data = load_integrated_empirical_data()
        
        for data_type, df in integrated_data.items():
            if df is not None:
                record_count = len(df)
                validation_results['data_quality_metrics'][data_type] = {
                    'records': record_count,
                    'status': 'integrated'
                }
                
                # Check for surveillance data indicators
                surveillance_indicators = 0
                if 'notes' in df.columns:
                    surveillance_records = df['notes'].str.contains('surveillance', na=False).sum()
                    surveillance_indicators = surveillance_records
                
                print(f"✅ {data_type:20} | {record_count:7,} total | {surveillance_indicators:6,} surveillance")
            else:
                validation_results['data_quality_metrics'][data_type] = {
                    'records': 0,
                    'status': 'not_available'
                }
                print(f"❌ {data_type:20} | No data available")
        
        print(f"\n🚀 Enhanced loader integration: SUCCESSFUL")
        validation_results['integration_status']['enhanced_loader'] = 'success'
        
    except Exception as e:
        print(f"❌ Enhanced loader integration failed: {e}")
        validation_results['integration_status']['enhanced_loader'] = f'error: {e}'
    
    # 3. Validate Plot Generation
    print("\n3. Plot Generation Validation")
    print("-" * 40)
    
    plot_categories = [
        'drug_usage_ddd_per_1000_per_day',
        'mean_mic_by_drug_per_bacteria', 
        'drug_failure_rate_by_bacteria_region',
        'mean_any_r_by_drug_for_each_bacteria'
    ]
    
    output_dir = Path('output_graphs')
    
    for category in plot_categories:
        category_dir = output_dir / category
        if category_dir.exists():
            plot_files = list(category_dir.glob('*.png'))
            validation_results['plot_generation_status'][category] = {
                'status': 'success',
                'plot_count': len(plot_files)
            }
            print(f"✅ {category:35} | {len(plot_files):4} plots generated")
        else:
            validation_results['plot_generation_status'][category] = {
                'status': 'missing'
            }
            print(f"❌ {category:35} | Directory not found")
    
    # 4. Coverage Assessment vs. Original System
    print("\n4. Coverage Assessment")
    print("-" * 40)
    
    # Compare with original ECDC data
    original_ecdc = Path('data/ecdc_resistance_2023.csv')
    if original_ecdc.exists():
        original_df = pd.read_csv(original_ecdc)
        original_records = len(original_df)
        
        enhanced_ecdc = Path('data/ecdc/ears_net_surveillance.csv')
        if enhanced_ecdc.exists():
            enhanced_df = pd.read_csv(enhanced_ecdc)
            enhanced_records = len(enhanced_df)
            
            improvement = (enhanced_records / original_records) if original_records > 0 else 0
            validation_results['coverage_assessment']['ecdc_improvement'] = {
                'original_records': original_records,
                'enhanced_records': enhanced_records,
                'improvement_factor': improvement
            }
            
            print(f"📈 ECDC Data Improvement:")
            print(f"   Original: {original_records:,} records")
            print(f"   Enhanced: {enhanced_records:,} records")
            print(f"   Improvement: {improvement:.1f}x")
    
    # Check empirical data summary
    summary_file = Path('data/empirical_data_summary.json')
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            summary = json.load(f)
            total_acquired = summary.get('total_records', 0)
            
            validation_results['coverage_assessment']['total_acquired_records'] = total_acquired
            print(f"\n🎯 Total empirical records acquired: {total_acquired:,}")
            print(f"📊 Sources successfully integrated: {len(summary['sources_acquired'])}")
    
    # 5. Final Assessment
    print("\n5. Final Assessment")
    print("-" * 40)
    
    successful_sources = len([s for s, status in validation_results['integration_status'].items() 
                             if isinstance(status, dict) and status.get('status') == 'success'])
    
    total_integrated_records = sum(
        metrics.get('records', 0) for metrics in validation_results['data_quality_metrics'].values()
        if isinstance(metrics, dict)
    )
    
    print(f"✅ Data sources successfully integrated: {successful_sources}/4")
    print(f"✅ Total integrated empirical records: {total_integrated_records:,}")
    print(f"✅ Enhanced empirical loader: {'Working' if validation_results['integration_status'].get('enhanced_loader') == 'success' else 'Failed'}")
    print(f"✅ Plot generation with empirical overlays: {'Working' if any(validation_results['plot_generation_status'].values()) else 'Failed'}")
    
    # Overall assessment
    if successful_sources >= 3 and total_integrated_records > 100000:
        overall_status = "EXCELLENT"
        status_icon = "🏆"
    elif successful_sources >= 2 and total_integrated_records > 50000:
        overall_status = "GOOD" 
        status_icon = "✅"
    elif successful_sources >= 1:
        overall_status = "PARTIAL"
        status_icon = "⚠️"
    else:
        overall_status = "FAILED"
        status_icon = "❌"
    
    validation_results['overall_assessment'] = {
        'status': overall_status,
        'successful_sources': successful_sources,
        'total_records': total_integrated_records,
        'timestamp': datetime.now().isoformat()
    }
    
    print(f"\n{status_icon} Overall Integration Status: {overall_status}")
    
    # Save validation results
    validation_file = Path('empirical_integration_validation.json')
    with open(validation_file, 'w') as f:
        json.dump(validation_results, f, indent=2)
    
    print(f"\n📋 Validation results saved to: {validation_file}")
    
    return validation_results

if __name__ == "__main__":
    results = validate_empirical_integration()
    
    # Print summary
    print("\n" + "="*60)
    print("🎉 EMPIRICAL DATA INTEGRATION VALIDATION COMPLETE")
    print("="*60)