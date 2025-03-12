from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, WheelZoomTool
from bokeh.embed import file_html
from bokeh.resources import CDN

def superimpose_plots(dict1, dict2, title="B-Y-Ion Mass Spectrum"):
    all_keys = sorted(set(dict1.keys()).union(dict2.keys()))
    dict1_vals = [dict1.get(k, 0) for k in all_keys]
    dict2_vals = [dict2.get(k, 0) for k in all_keys]


    source = ColumnDataSource(data=dict(
        x_values=[x for x in all_keys],
        reds=dict1_vals,
        blues=dict2_vals,
        labels_x = [k if dict1.get(k, 0) > 0 else " " for k in all_keys],
        labels_y = [k if dict2.get(k, 0) > 0 else " " for k in all_keys]
    ))

    p = figure(x_range=(min(all_keys)-20, max(all_keys)+20), title=title, toolbar_location=None,
               tools="xpan,reset,save", x_axis_label = "Mass M/z", y_axis_label = "Intensity" )

    wheel_zoom_tool = WheelZoomTool(dimensions='width')
    p.add_tools(wheel_zoom_tool)
    p.toolbar.active_scroll = wheel_zoom_tool

    hover = HoverTool(tooltips=[
        ("Mass", "@x_values"),
        ("B-Ions", "@reds"),
        ("Y-Ions", "@blues")
    ])
    p.add_tools(hover)

    p.vbar(x='x_values', top='reds', width=0.4, source=source, legend_label="B-Ions", color="red", muted_color="red",
           muted_alpha=0.2)
    p.vbar(x='x_values', top='blues', width=0.4, source=source, legend_label="Y-Ions", color="blue", muted_color="blue",
           muted_alpha=0.2)

    labels_x = LabelSet(x='x_values', y='reds', text='labels_x', level='glyph',
                      x_offset=10, y_offset=5, source=source,
                      text_font_size="6pt", text_color="red")
    p.add_layout(labels_x)

    labels_y = LabelSet(x='x_values', y='blues', text='labels_y', level='glyph',
                      x_offset=10, y_offset=0, source=source,
                      text_font_size="6pt", text_color="blue")
    p.add_layout(labels_y)

    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    p.legend.location = "top_right"
    p.legend.click_policy = "mute"

    html = file_html(p, CDN)

    return html


if __name__ == '__main__':
    ms = {802.4007260000001: 0.6204353503100348, 803.400921: 0.24436153875197233, 804.4029763333333: 0.08038456698870841, 805.399877: 0.010007372238203072}
    mz = {702.4007260000001: 0.6204353503100348, 703.400921: 0.24436153875197233, 704.4029763333333: 0.08038456698870841, 705.399877: 0.010007372238203072}

    dd = superimpose_plots(ms, mz)