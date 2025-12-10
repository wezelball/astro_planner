import plotly.graph_objects as go
import numpy as np

def plot_sky_polar(df, horizon):
    """
    Polar sky plot with object positions and shaded local horizon.
    Azimuth = theta, Altitude = r
    """

    import plotly.graph_objects as go

    fig = go.Figure()

    # Bigger default figure and looser margins so outer labels are not clipped
    fig.update_layout(
        width=1000,
        height=800,
        margin=dict(l=80, r=80, t=80, b=80),
    )

    # ----------------------------
    # Horizon processing (FIX GAP + SHADE)
    # ----------------------------
    if horizon:
        hz_az = [pt[0] for pt in horizon]
        hz_alt = [pt[1] for pt in horizon]

        # Ensure wrap-around at 0° and 360°
        if hz_az[0] != 0.0:
            hz_az.insert(0, 0.0)
            hz_alt.insert(0, hz_alt[0])

        if hz_az[-1] != 360.0:
            hz_az.append(360.0)
            hz_alt.append(hz_alt[0])

        # Horizon line
        fig.add_trace(go.Scatterpolar(
            theta=hz_az,
            r=hz_alt,
            mode="lines",
            line=dict(width=2),
            name="Horizon",
        ))

        # Shaded horizon (everything below horizon)
        fig.add_trace(go.Scatterpolar(
            theta=hz_az,
            r=hz_alt,
            fill="toself",
            fillcolor="rgba(100,100,100,0.4)",  # semi-transparent gray
            line=dict(width=0),
            name="Blocked Sky",
            showlegend=True
        ))

    # ----------------------------------------
    # Add cardinal & intercardinal direction labels (outside plotted range)
    # ----------------------------------------
    # Choose a radius for the labels that's outside the max altitude (90)
    label_radius = 100  # put labels outside the visible 0-90 altitude range
    label_font = dict(size=16, color="white", family="Arial")

    cardinal_labels = {
        0:  "N",
        45: "NE",
        90: "E",
        135:"SE",
        180:"S",
        225:"SW",
        270:"W",
        315:"NW",
    }

    for az, label in cardinal_labels.items():
        fig.add_trace(
            go.Scatterpolar(
                r=[label_radius],           # outside the 0-90 altitude range
                theta=[az],
                mode="text",
                text=[label],
                textposition="middle center",
                textfont=label_font,
                showlegend=False,
                hoverinfo="skip"
            )
        )

    # ----------------------------
    # Plot objects
    # ----------------------------
    if df is not None and not df.empty:
        hover_text = [
            f"<b>{name}</b><br>Az: {az:.1f}°<br>Alt: {alt:.1f}°"
            for name, az, alt in zip(df["name"], df["az_deg"], df["alt_deg"])
        ]

        fig.add_trace(go.Scatterpolar(
            theta=df["az_deg"],
            r=df["alt_deg"],
            mode="markers+text",
            text=df["name"],
            textposition="top center",
            textfont=dict(color="red", size=11),
            marker=dict(size=8),
            hovertext=hover_text,      # custom tooltip
            hoverinfo="text",          # hide r/theta default display
            name="Objects"
        ))

    # ----------------------------
    # Layout
    # ----------------------------
    # Set radial axis max to include label_radius (give a bit of margin)
    radial_max = max(95, label_radius + 5)

    fig.update_layout(
        title="Sky View with Local Horizon",
        polar=dict(
            radialaxis=dict(
                range=[0, radial_max],
                tickangle=90,
                title="Altitude (deg)",
                ticks="",           # hides tick marks
                showticklabels=False
            ),
            angularaxis=dict(
                direction="clockwise",
                rotation=90,
                tickmode="array",
                tickvals=[0,45,90,135,180,225,270,315],   # optional axis ticks
                ticktext=["N","NE","E","SE","S","SW","W","NW"]
            )
        ),
        showlegend=True,
        margin=dict(l=80, r=80, t=80, b=80)
    )

    return fig
