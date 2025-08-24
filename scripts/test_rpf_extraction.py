#!/usr/bin/env python3
"""Comprehensive testing of RPF extraction system.

Tests the RPF extraction system on all generated test cases and provides
detailed performance analysis.
"""

import sys
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from getRPF.core.processors.rpf_extractor import RPFExtractor

def test_single_case(
    extractor: RPFExtractor, 
    test_case: str, 
    metadata: Dict
) -> Dict:
    """Test RPF extraction on a single test case."""
    
    input_file = Path(metadata['file_path'])
    output_file = Path(f"comprehensive_test_data/{test_case}_extracted.fastq")
    format_type = 'collapsed' if metadata.get('format') == 'collapsed' else 'fastq'
    
    print(f"🧪 Testing {test_case}...")
    print(f"   Description: {metadata['description']}")
    
    start_time = time.time()
    
    try:
        results = extractor.extract_rpfs(
            input_file=input_file,
            output_file=output_file,
            format=format_type,
            max_reads=None  # Process all reads
        )
        
        processing_time = time.time() - start_time
        
        # Calculate metrics
        actual_success_rate = results.extracted_rpfs / results.input_reads if results.input_reads > 0 else 0
        expected_success_rate = metadata['expected_success_rate']
        
        # Check if method matches expectation
        method_match = results.extraction_method == metadata.get('expected_method', 'unknown')
        
        # Architecture match (if expected)
        arch_match = True
        if 'expected_architecture' in metadata:
            arch_match = results.architecture_match == metadata['expected_architecture']
        
        test_result = {
            'test_case': test_case,
            'status': 'PASS' if actual_success_rate >= expected_success_rate * 0.9 else 'FAIL',
            'input_reads': results.input_reads,
            'extracted_rpfs': results.extracted_rpfs,
            'actual_success_rate': actual_success_rate,
            'expected_success_rate': expected_success_rate,
            'extraction_method': results.extraction_method,
            'expected_method': metadata.get('expected_method', 'N/A'),
            'method_match': method_match,
            'architecture_match': results.architecture_match,
            'expected_architecture': metadata.get('expected_architecture', 'N/A'),
            'arch_match': arch_match,
            'processing_time': processing_time,
            'quality_metrics': results.quality_metrics,
            'error': None
        }
        
        # Print results
        success_pct = actual_success_rate * 100
        expected_pct = expected_success_rate * 100
        status_emoji = "✅" if test_result['status'] == 'PASS' else "❌"
        
        print(f"   {status_emoji} {test_result['status']}: {success_pct:.1f}% success (expected {expected_pct:.0f}%)")
        print(f"   Method: {results.extraction_method} {'✓' if method_match else '✗'}")
        print(f"   Extracted: {results.extracted_rpfs}/{results.input_reads} reads")
        print(f"   Time: {processing_time:.2f}s")
        
    except Exception as e:
        test_result = {
            'test_case': test_case,
            'status': 'ERROR',
            'error': str(e),
            'processing_time': time.time() - start_time
        }
        print(f"   ❌ ERROR: {e}")
    
    return test_result

def analyze_results(results: List[Dict]) -> None:
    """Analyze and print overall test results."""
    
    print("\n" + "="*80)
    print("📊 COMPREHENSIVE TEST RESULTS ANALYSIS")
    print("="*80)
    
    # Overall statistics
    total_tests = len(results)
    passed_tests = sum(1 for r in results if r['status'] == 'PASS')
    failed_tests = sum(1 for r in results if r['status'] == 'FAIL')
    error_tests = sum(1 for r in results if r['status'] == 'ERROR')
    
    print(f"\n🎯 Overall Performance:")
    print(f"   Total tests: {total_tests}")
    print(f"   Passed: {passed_tests} ({passed_tests/total_tests*100:.1f}%)")
    print(f"   Failed: {failed_tests} ({failed_tests/total_tests*100:.1f}%)")
    print(f"   Errors: {error_tests} ({error_tests/total_tests*100:.1f}%)")
    
    # Success rate analysis
    valid_results = [r for r in results if r['status'] in ['PASS', 'FAIL']]
    if valid_results:
        total_reads = sum(r['input_reads'] for r in valid_results)
        total_extracted = sum(r['extracted_rpfs'] for r in valid_results)
        overall_success_rate = total_extracted / total_reads if total_reads > 0 else 0
        
        print(f"\n📈 Extraction Performance:")
        print(f"   Total reads processed: {total_reads:,}")
        print(f"   Total RPFs extracted: {total_extracted:,}")
        print(f"   Overall success rate: {overall_success_rate*100:.2f}%")
    
    # Method analysis
    method_counts = {}
    for r in valid_results:
        method = r.get('extraction_method', 'unknown')
        method_counts[method] = method_counts.get(method, 0) + 1
    
    print(f"\n🔧 Extraction Methods Used:")
    for method, count in method_counts.items():
        print(f"   {method}: {count} tests ({count/len(valid_results)*100:.1f}%)")
    
    # Performance by test case
    print(f"\n📋 Detailed Results by Test Case:")
    print(f"{'Test Case':<25} {'Status':<6} {'Success Rate':<12} {'Method':<15} {'Time':<8}")
    print("-" * 80)
    
    for r in results:
        if r['status'] == 'ERROR':
            print(f"{r['test_case']:<25} {'ERROR':<6} {'N/A':<12} {'N/A':<15} {r['processing_time']:.2f}s")
        else:
            success_rate = f"{r['actual_success_rate']*100:.1f}%"
            method = r['extraction_method'][:14]
            print(f"{r['test_case']:<25} {r['status']:<6} {success_rate:<12} {method:<15} {r['processing_time']:.2f}s")
    
    # Identify problem areas
    failed_results = [r for r in results if r['status'] == 'FAIL']
    if failed_results:
        print(f"\n⚠️  Failed Test Cases (need improvement):")
        for r in failed_results:
            expected = r['expected_success_rate'] * 100
            actual = r['actual_success_rate'] * 100
            gap = expected - actual
            print(f"   {r['test_case']}: {actual:.1f}% vs {expected:.0f}% expected (-{gap:.1f}%)")
    
    # Success stories
    excellent_results = [r for r in valid_results if r['actual_success_rate'] >= 0.95]
    if excellent_results:
        print(f"\n🏆 Excellent Performance (≥95% success):")
        for r in excellent_results:
            print(f"   {r['test_case']}: {r['actual_success_rate']*100:.1f}% success with {r['extraction_method']}")

def main():
    """Run comprehensive RPF extraction tests."""
    
    print("🚀 COMPREHENSIVE RPF EXTRACTION TESTING")
    print("="*50)
    
    # Load test metadata
    metadata_file = Path("comprehensive_test_data/test_cases_metadata.json")
    if not metadata_file.exists():
        print("❌ ERROR: Test cases not found. Run generate_test_cases.py first.")
        return
    
    with open(metadata_file) as f:
        all_metadata = json.load(f)
    
    print(f"Found {len(all_metadata)} test cases to evaluate...")
    
    # Initialize extractor
    extractor = RPFExtractor()
    
    # Run tests
    all_results = []
    start_time = time.time()
    
    for test_case, metadata in all_metadata.items():
        result = test_single_case(extractor, test_case, metadata)
        all_results.append(result)
        print()  # Empty line between tests
    
    total_time = time.time() - start_time
    
    # Analyze results
    analyze_results(all_results)
    
    print(f"\n⏱️  Total testing time: {total_time:.1f}s")
    
    # Save detailed results
    results_file = Path("comprehensive_test_data/test_results.json")
    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': time.time(),
            'total_tests': len(all_results),
            'total_time': total_time,
            'results': all_results
        }, f, indent=2)
    
    print(f"💾 Detailed results saved to {results_file}")
    
    # Final verdict
    passed_tests = sum(1 for r in all_results if r['status'] == 'PASS')
    total_tests = len(all_results)
    
    if passed_tests == total_tests:
        print(f"\n🎉 ALL TESTS PASSED! RPF extraction system is robust and ready for production.")
    elif passed_tests / total_tests >= 0.8:
        print(f"\n✅ GOOD PERFORMANCE: {passed_tests}/{total_tests} tests passed. System is reliable with room for minor improvements.")
    else:
        print(f"\n⚠️  NEEDS IMPROVEMENT: {passed_tests}/{total_tests} tests passed. Consider refining algorithms for edge cases.")

if __name__ == "__main__":
    main()