try:
    from Bio import SeqIO
    import pandas as pd
    import requests
    import plotly.graph_objects as go
    print("All libraries installed successfully!")
except ImportError as e:
    print(f"x Error: {e}")

