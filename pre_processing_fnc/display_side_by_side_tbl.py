from IPython.display import display, HTML
import pandas as pd

def display_side_by_side(dfs, titles=[]):
    """
    Display DataFrames side by side in Jupyter Notebook.
    
    Parameters:
    dfs (list): List of DataFrames or tuples to display.
    titles (list): List of titles for the DataFrames (optional).
    """
    html_str = ''
    
    for i, df in enumerate(dfs):
        # Convert tuples to DataFrames if necessary
        if isinstance(df, tuple):
            df = pd.DataFrame(df)
        elif not isinstance(df, pd.DataFrame):
            raise TypeError(f"Expected pd.DataFrame or tuple, but got {type(df)} at index {i}")
        
        title = f'<h3>{titles[i]}</h3>' if i < len(titles) else ''
        html_str += f'<td>{title}{df.to_html()}</td>'
    
    display(HTML(f"""
    <table>
        <tr>{html_str}</tr>
    </table>
    """))