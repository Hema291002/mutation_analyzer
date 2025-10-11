#!/usr/bin/env python3
"""
Mutation Analysis Visualization Generator
Creates interactive HTML visualizations using Plotly
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import os

def create_sequence_viewer(original_seq, mutated_seq, position, gene_name, mutation_notation):
    """Create interactive sequence comparison visualization"""
    
    # Show window around mutation
    window = 30
    start = max(0, position - window)
    end = min(len(original_seq), position + window + 1)
    
    positions = list(range(start, end))
    orig_bases = list(original_seq[start:end])
    mut_bases = list(mutated_seq[start:end])
    
    # Color coding
    colors = ['red' if pos == position else 'lightblue' for pos in positions]
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Original Sequence', 'Mutated Sequence'),
        vertical_spacing=0.15
    )
    
    # Original sequence
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers+text',
            marker=dict(size=20, color='lightblue', line=dict(width=2, color='navy')),
            text=orig_bases,
            textposition='middle center',
            textfont=dict(size=12, color='black', family='Courier'),
            name='Original',
            hovertemplate='<b>Position %{x}</b><br>Base: %{text}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Mutated sequence
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers+text',
            marker=dict(size=20, color=colors, line=dict(width=2, color='darkred')),
            text=mut_bases,
            textposition='middle center',
            textfont=dict(size=12, color='black', family='Courier'),
            name='Mutated',
            hovertemplate='<b>Position %{x}</b><br>Base: %{text}<extra></extra>'
        ),
        row=2, col=1
    )
    
    fig.update_layout(
        title=dict(
            text=f'<b>Sequence Visualization: {gene_name} - {mutation_notation}</b>',
            font=dict(size=20)
        ),
        showlegend=False,
        height=500,
        hovermode='closest',
        plot_bgcolor='#f8f9fa'
    )
    
    fig.update_xaxes(title_text="Genomic Position", row=2, col=1)
    fig.update_yaxes(showticklabels=False)
    
    return fig


def create_codon_visualization(classification):
    """Create codon change comparison"""
    
    orig_codon = classification['original_codon']
    mut_codon = classification['mutated_codon']
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            f'Original: {orig_codon} → {classification["ref_aa"]}',
            f'Mutated: {mut_codon} → {classification["alt_aa"]}'
        ),
        horizontal_spacing=0.15
    )
    
    positions = [0, 1, 2]
    
    # Original codon
    colors_orig = ['#4CAF50' if orig_codon[i] == mut_codon[i] else '#FF9800' for i in range(3)]
    
    fig.add_trace(
        go.Bar(
            x=positions,
            y=[1, 1, 1],
            marker=dict(color=colors_orig, line=dict(width=2, color='black')),
            text=list(orig_codon),
            textposition='inside',
            textfont=dict(size=32, color='white', family='Courier'),
            name='Original',
            hovertemplate='Position %{x}<br>Base: %{text}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Mutated codon
    colors_mut = ['#4CAF50' if orig_codon[i] == mut_codon[i] else '#F44336' for i in range(3)]
    
    fig.add_trace(
        go.Bar(
            x=positions,
            y=[1, 1, 1],
            marker=dict(color=colors_mut, line=dict(width=2, color='black')),
            text=list(mut_codon),
            textposition='inside',
            textfont=dict(size=32, color='white', family='Courier'),
            name='Mutated',
            hovertemplate='Position %{x}<br>Base: %{text}<extra></extra>'
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        title=dict(
            text=f'<b>Codon Change Analysis: {classification["notation"]}</b>',
            font=dict(size=20)
        ),
        showlegend=False,
        height=400,
        plot_bgcolor='#f8f9fa'
    )
    
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    
    return fig


def create_pathogenicity_gauge(pathogenicity):
    """Create gauge chart for pathogenicity score"""
    
    score = pathogenicity['score']
    risk = pathogenicity['risk']
    
    # Color based on risk
    if risk == 'HIGH':
        bar_color = '#F44336'
    elif risk == 'MODERATE':
        bar_color = '#FF9800'
    else:
        bar_color = '#4CAF50'
    
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=score,
        domain={'x': [0, 1], 'y': [0, 1]},
        title={'text': f"<b>Pathogenicity Score</b><br><span style='font-size:0.8em'>{risk} RISK</span>", 
               'font': {'size': 24}},
        delta={'reference': 2.5, 'increasing': {'color': "red"}},
        gauge={
            'axis': {'range': [None, 5], 'tickwidth': 2, 'tickcolor': "black"},
            'bar': {'color': bar_color, 'thickness': 0.75},
            'bgcolor': "white",
            'borderwidth': 3,
            'bordercolor': "gray",
            'steps': [
                {'range': [0, 1.5], 'color': '#C8E6C9'},
                {'range': [1.5, 3.5], 'color': '#FFE0B2'},
                {'range': [3.5, 5], 'color': '#FFCDD2'}
            ],
            'threshold': {
                'line': {'color': "black", 'width': 4},
                'thickness': 0.75,
                'value': score
            }
        }
    ))
    
    fig.update_layout(
        height=400,
        margin=dict(l=20, r=20, t=80, b=20),
        font={'family': 'Arial'}
    )
    
    return fig


def create_mutation_type_chart(classification):
    """Create mutation type indicator"""
    
    mutation_types = {
        'Missense': {'color': '#FF9800', 'impact': 'Moderate', 'icon': '⚠️'},
        'Nonsense': {'color': '#F44336', 'impact': 'High', 'icon': '🛑'},
        'Silent/Synonymous': {'color': '#4CAF50', 'impact': 'Low', 'icon': '✓'},
        'Stop-loss': {'color': '#E91E63', 'impact': 'High', 'icon': '⚠️'}
    }
    
    current_type = classification['type']
    type_info = mutation_types.get(current_type, {'color': '#9E9E9E', 'impact': 'Unknown', 'icon': '?'})
    
    # Create a simple annotation-based display instead
    fig = go.Figure()
    
    fig.add_annotation(
        text=f"<b>{type_info['icon']}<br>{current_type}</b><br><br>" +
             f"<span style='font-size:18px'>Impact Level: {type_info['impact']}</span>",
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        showarrow=False,
        font=dict(size=28, color='white', family='Arial'),
        align='center'
    )
    
    fig.update_layout(
        height=300,
        paper_bgcolor=type_info['color'],
        xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
        yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
        margin=dict(l=20, r=20, t=20, b=20),
        plot_bgcolor=type_info['color']
    )
    
    return fig


def create_pathway_sunburst(pathways, gene_name):
    """Create sunburst chart for pathway visualization"""
    
    if not pathways:
        return None
    
    labels = [gene_name]
    parents = ['']
    values = [len(pathways)]
    colors = ['#1f77b4']
    
    color_palette = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
    
    for i, (pathway, _) in enumerate(pathways):
        labels.append(pathway)
        parents.append(gene_name)
        values.append(1)
        colors.append(color_palette[i % len(color_palette)])
    
    fig = go.Figure(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        marker=dict(colors=colors, line=dict(width=2, color='white')),
        hovertemplate='<b>%{label}</b><br>%{value} pathway(s)<extra></extra>',
        textfont=dict(size=14)
    ))
    
    fig.update_layout(
        title=dict(
            text=f'<b>Biological Pathways: {gene_name}</b>',
            font=dict(size=20)
        ),
        height=600,
        margin=dict(t=100, l=0, r=0, b=0)
    )
    
    return fig


def create_disease_chart(cosmic, omim, gene_name):
    """Create disease association bar chart"""
    
    diseases = []
    sources = []
    colors_map = []
    
    if omim and omim.get('diseases'):
        for disease in omim['diseases']:
            diseases.append(disease)
            sources.append('OMIM (Genetic)')
            colors_map.append('#E74C3C')
    
    if cosmic and cosmic.get('cancer_types'):
        for cancer in cosmic['cancer_types'][:4]:
            diseases.append(cancer)
            sources.append('COSMIC (Cancer)')
            colors_map.append('#3498DB')
    
    if not diseases:
        return None
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        y=diseases,
        x=[1] * len(diseases),
        orientation='h',
        marker=dict(color=colors_map, line=dict(width=2, color='black')),
        text=sources,
        textposition='inside',
        textfont=dict(color='white', size=12),
        hovertemplate='<b>%{y}</b><br>Source: %{text}<extra></extra>'
    ))
    
    fig.update_layout(
        title=dict(
            text=f'<b>Disease Associations: {gene_name}</b>',
            font=dict(size=20)
        ),
        xaxis=dict(showticklabels=False, title=''),
        yaxis=dict(title='Disease/Condition', tickfont=dict(size=12)),
        height=400,
        plot_bgcolor='#f8f9fa',
        showlegend=False
    )
    
    return fig


def create_comprehensive_dashboard(gene_name, sequence, mutated_seq, position, 
                                   classification, pathogenicity, cosmic, omim, pathways):
    """Create complete HTML dashboard with all visualizations"""
    
    # Generate all figures
    seq_fig = create_sequence_viewer(sequence, mutated_seq, position, gene_name, classification['notation'])
    codon_fig = create_codon_visualization(classification)
    gauge_fig = create_pathogenicity_gauge(pathogenicity)
    type_fig = create_mutation_type_chart(classification)
    pathway_fig = create_pathway_sunburst(pathways, gene_name) if pathways else None
    disease_fig = create_disease_chart(cosmic, omim, gene_name)
    
    # Create HTML dashboard
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Mutation Analysis Dashboard - {gene_name}</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                padding: 20px;
                min-height: 100vh;
            }}
            
            .container {{
                max-width: 1400px;
                margin: 0 auto;
            }}
            
            .header {{
                background: white;
                padding: 40px;
                border-radius: 15px;
                margin-bottom: 30px;
                box-shadow: 0 10px 30px rgba(0,0,0,0.3);
                text-align: center;
            }}
            
            .header h1 {{
                font-size: 2.5em;
                color: #667eea;
                margin-bottom: 10px;
            }}
            
            .header .subtitle {{
                font-size: 1.3em;
                color: #666;
            }}
            
            .summary-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px;
                margin-bottom: 30px;
            }}
            
            .card {{
                background: white;
                padding: 25px;
                border-radius: 12px;
                box-shadow: 0 5px 15px rgba(0,0,0,0.2);
                transition: transform 0.3s ease;
            }}
            
            .card:hover {{
                transform: translateY(-5px);
            }}
            
            .card h3 {{
                color: #667eea;
                font-size: 1.1em;
                margin-bottom: 10px;
                text-transform: uppercase;
                letter-spacing: 1px;
            }}
            
            .card .value {{
                font-size: 2em;
                font-weight: bold;
                color: #333;
            }}
            
            .card .label {{
                font-size: 0.9em;
                color: #999;
                margin-top: 5px;
            }}
            
            .risk-high {{ color: #F44336; }}
            .risk-moderate {{ color: #FF9800; }}
            .risk-low {{ color: #4CAF50; }}
            
            .viz-section {{
                background: white;
                padding: 30px;
                border-radius: 12px;
                margin-bottom: 20px;
                box-shadow: 0 5px 15px rgba(0,0,0,0.2);
            }}
            
            .viz-section h2 {{
                color: #667eea;
                margin-bottom: 20px;
                padding-bottom: 10px;
                border-bottom: 3px solid #667eea;
            }}
            
            .two-col {{
                display: grid;
                grid-template-columns: 1fr 1fr;
                gap: 20px;
            }}
            
            @media (max-width: 768px) {{
                .two-col {{
                    grid-template-columns: 1fr;
                }}
            }}
            
            .footer {{
                text-align: center;
                color: white;
                padding: 20px;
                margin-top: 30px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 Mutation Analysis Dashboard</h1>
                <p class="subtitle">Gene: <strong>{gene_name}</strong> | Variant: <strong>{classification['notation']}</strong></p>
            </div>
            
            <div class="summary-grid">
                <div class="card">
                    <h3>Mutation Type</h3>
                    <div class="value">{classification['type']}</div>
                    <div class="label">Classification</div>
                </div>
                
                <div class="card">
                    <h3>Pathogenicity</h3>
                    <div class="value risk-{pathogenicity['risk'].lower()}">{pathogenicity['prediction']}</div>
                    <div class="label">Clinical Significance</div>
                </div>
                
                <div class="card">
                    <h3>Risk Level</h3>
                    <div class="value risk-{pathogenicity['risk'].lower()}">{pathogenicity['risk']}</div>
                    <div class="label">Assessment</div>
                </div>
                
                <div class="card">
                    <h3>Score</h3>
                    <div class="value">{pathogenicity['score']}/5</div>
                    <div class="label">Pathogenicity Index</div>
                </div>
            </div>
            
            <div class="viz-section">
                <h2>📊 Sequence Visualization</h2>
                <div id="seq-plot"></div>
            </div>
            
            <div class="two-col">
                <div class="viz-section">
                    <h2>🧬 Codon Analysis</h2>
                    <div id="codon-plot"></div>
                </div>
                
                <div class="viz-section">
                    <h2>📈 Mutation Type</h2>
                    <div id="type-plot"></div>
                </div>
            </div>
            
            <div class="viz-section">
                <h2>⚕️ Pathogenicity Assessment</h2>
                <div id="gauge-plot"></div>
            </div>
            
            {f'''
            <div class="viz-section">
                <h2>🔬 Biological Pathways</h2>
                <div id="pathway-plot"></div>
            </div>
            ''' if pathway_fig else ''}
            
            {f'''
            <div class="viz-section">
                <h2>🏥 Disease Associations</h2>
                <div id="disease-plot"></div>
            </div>
            ''' if disease_fig else ''}
            
            <div class="viz-section">
                <h2>📋 Detailed Information</h2>
                <table style="width: 100%; border-collapse: collapse;">
                    <tr style="background: #f5f5f5;">
                        <td style="padding: 15px; border: 1px solid #ddd;"><strong>DNA Change</strong></td>
                        <td style="padding: 15px; border: 1px solid #ddd;">{classification['original_codon']} → {classification['mutated_codon']}</td>
                    </tr>
                    <tr>
                        <td style="padding: 15px; border: 1px solid #ddd;"><strong>Protein Change</strong></td>
                        <td style="padding: 15px; border: 1px solid #ddd;">{classification['ref_aa']} → {classification['alt_aa']}</td>
                    </tr>
                    <tr style="background: #f5f5f5;">
                        <td style="padding: 15px; border: 1px solid #ddd;"><strong>Position</strong></td>
                        <td style="padding: 15px; border: 1px solid #ddd;">Codon {classification['codon_number']} (Nucleotide {position + 1})</td>
                    </tr>
                    <tr>
                        <td style="padding: 15px; border: 1px solid #ddd;"><strong>Recommendation</strong></td>
                        <td style="padding: 15px; border: 1px solid #ddd;">{pathogenicity['recommendation']}</td>
                    </tr>
                </table>
            </div>
            
            <div class="footer">
                <p>Generated by Enhanced Mutation Analyzer | {gene_name} Analysis</p>
                <p style="font-size: 0.9em; margin-top: 10px;">For research and educational purposes only</p>
            </div>
        </div>
        
        <script>
            var seqData = {seq_fig.to_json()};
            Plotly.newPlot('seq-plot', seqData.data, seqData.layout, {{responsive: true}});
            
            var codonData = {codon_fig.to_json()};
            Plotly.newPlot('codon-plot', codonData.data, codonData.layout, {{responsive: true}});
            
            var gaugeData = {gauge_fig.to_json()};
            Plotly.newPlot('gauge-plot', gaugeData.data, gaugeData.layout, {{responsive: true}});
            
            var typeData = {type_fig.to_json()};
            Plotly.newPlot('type-plot', typeData.data, typeData.layout, {{responsive: true}});
            
            {f'''
            var pathwayData = {pathway_fig.to_json()};
            Plotly.newPlot('pathway-plot', pathwayData.data, pathwayData.layout, {{responsive: true}});
            ''' if pathway_fig else ''}
            
            {f'''
            var diseaseData = {disease_fig.to_json()};
            Plotly.newPlot('disease-plot', diseaseData.data, diseaseData.layout, {{responsive: true}});
            ''' if disease_fig else ''}
        </script>
    </body>
    </html>
    """
    
    # Save dashboard
    os.makedirs('output', exist_ok=True)
    filename = f"output/{gene_name}_{classification['notation']}_dashboard.html"
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n🎨 Interactive dashboard created!")
    print(f"   📁 Saved to: {filename}")
    print(f"\n🌐 Open this file in your web browser to view the interactive visualizations!")
    
    return filename


if __name__ == "__main__":
    print("This is a visualization module.")
    print("Use it with analyzer_enhanced.py or import the functions.")