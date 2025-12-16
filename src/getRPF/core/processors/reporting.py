"""Reporting engine for RPF extraction.

Generates:
1. Rich CLI reports with ASCII sparklines.
2. Lightweight, interactive HTML reports (FastQC-style).
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import asdict

from .signals import SignalStats
from .types import SegmentInfo, ReadArchitecture


class Reporter:
    """Generates human-readable interaction reports."""

    def generate_cli_report(self, 
                            architecture: ReadArchitecture, 
                            segments: List[SegmentInfo], 
                            stats: SignalStats) -> str:
        """Generate a Rich-text compatible string for CLI output."""
        report = []
        report.append(f"✓ Architecture Detect: [bold green]{architecture.protocol_name}[/bold green]")
        
        # Structure Map
        structure_str = ""
        for seg in segments:
            length = seg.end_pos - seg.start_pos
            name = seg.segment_type.upper()
            structure_str += f"[{name}:{length}]--"
            
        report.append(f"  Structure: {structure_str[:-2]}") # Remove last --
        
        # ASCII Sparkline for Entropy
        sparkline = self._ascii_sparkline(stats.entropy_5p[:50])
        report.append(f"  Entropy Signal (5'): {sparkline}")
        
        return "\n".join(report)

    def generate_html_report(self, 
                             output_path: Path,
                             architecture: ReadArchitecture,
                             segments: List[SegmentInfo],
                             stats: SignalStats,
                             trace_log: List[str]):
        """Generate a single-file interactive HTML report."""
        
        # Prepare data for JS embedding
        entropy_data = list(stats.entropy_5p)
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>getRPF Extraction Report</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; color: #333; }}
        .header {{ border-bottom: 2px solid #eaeaea; padding-bottom: 20px; margin-bottom: 30px; }}
        .card {{ background: #fff; border: 1px solid #eaeaea; border-radius: 8px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        .metric-plot {{ width: 100%; height: 300px; background: #fafafa; border: 1px dashed #ddd; position: relative; }}
        .bar {{ display: inline-block; width: 10px; background: #0070f3; position: absolute; bottom: 0; }}
        h2 {{ margin-top: 0; font-size: 1.2rem; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
        .log-box {{ background: #f5f5f5; padding: 10px; border-radius: 4px; font-family: monospace; max-height: 200px; overflow-y: auto; }}
        .structure-map {{ display: flex; gap: 5px; margin: 20px 0; }}
        .segment {{ padding: 10px; border-radius: 4px; text-align: center; color: white; font-weight: bold; min-width: 50px; }}
        .umi {{ background: #f5a623; }}
        .rpf {{ background: #4a90e2; flex-grow: 1; }}
        .adapter {{ background: #d0021b; }}
        .barcode {{ background: #bd10e0; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>getRPF Extraction Report</h1>
        <p>File: {output_path.name}</p>
        <p>Identified Protocol: <strong>{architecture.protocol_name}</strong></p>
    </div>

    <div class="card">
        <h2>Read Structure</h2>
        <div class="structure-map">
            {"".join(self._render_html_segments(segments))}
        </div>
    </div>

    <div class="card">
        <h2>Signal Metrics (5' Entropy)</h2>
        <div class="metric-plot" id="entropyPlot">
            <!-- Simple JS/CSS Bars for now -->
            {self._render_css_bars(entropy_data[:50])}
        </div>
        <p><small>Showing first 50 positions. High entropy = Random (UMI), Low = Consensus (Adapter)</small></p>
    </div>

    <div class="card">
        <h2>Decision Trace</h2>
        <div class="log-box">
            { "<br>".join(trace_log) }
        </div>
    </div>
</body>
</html>
        """
        
        with open(output_path, "w") as f:
            f.write(html_content)

    def _ascii_sparkline(self, data: List[float]) -> str:
        """Create a simple ASCII sparkline."""
        chars = "  ▂▄▆█"
        if not data:
            return ""
        min_val, max_val = min(data), max(data)
        if max_val == min_val:
            max_val += 1e-6
            
        spark = ""
        for x in data:
            idx = int((x - min_val) / (max_val - min_val) * (len(chars) - 1))
            spark += chars[idx]
        return spark

    def _render_html_segments(self, segments: List[SegmentInfo]) -> List[str]:
        """Render HTML divs for segments."""
        html_segs = []
        for seg in segments:
            type_class = seg.segment_type.lower()
            length = seg.end_pos - seg.start_pos
            html_segs.append(f'<div class="segment {type_class}">{seg.segment_type.upper()}<br>{length}nt</div>')
        return html_segs

    def _render_css_bars(self, data: List[float]) -> str:
        """Render simple CSS bars for the plot."""
        bars = ""
        max_val = max(data) if data else 1
        width_pct = 100 / len(data) if data else 0
        
        for i, val in enumerate(data):
            height_pct = (val / max_val) * 100
            left_pct = i * width_pct
            bars += f'<div class="bar" style="left: {left_pct}%; width: {width_pct}%; height: {height_pct}%;" title="Pos {i}: {val:.2f}"></div>'
        return bars
