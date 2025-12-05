import plotly.graph_objects as go
import numpy as np

def plot_sky_polar(df, horizon):
    """
    Polar sky plot with object positions and shaded local horizon.
    Azimuth = theta, Altitude = r
    """

    fig = go.Figure()

    fig.update_layout(
    width=950,     # was likely default ~700
    height=750,    # increases vertical resolution
    margin=dict(l=60, r=60, t=60, b=60)
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

    # ----------------------------
    # Plot objects
    # ----------------------------
    if not df.empty:
        hover_text = [
            #f"<b>{name}</b><br>Alt: {az:.1f}°<br>Az: {alt:.1f}°"
            f"<b>{name}</b><br>Az: {az:.1f}°<br>Alt: {alt:.1f}°"
            #for name, alt, az in zip(df["name"], df["alt_deg"], df["az_deg"])
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
            hovertext=hover_text,      # ✅ CUSTOM TOOLTIP
            hoverinfo="text",          # ✅ HIDE r/theta
            name="Objects"
        ))

    # ----------------------------
    # Layout
    # ----------------------------
    fig.update_layout(
        title="Sky View with Local Horizon",
        polar=dict(
            radialaxis=dict(
                range=[0, 90],
                tickangle=90,
                title="Altitude (deg)"
            ),
            angularaxis=dict(
                direction="clockwise",
                rotation=90
            )
        ),
        showlegend=True,
        margin=dict(l=40, r=40, t=60, b=40)
    )

    return fig
