import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

ROLE_COLORS = {
    "Hedef ARG (Direnç Geni)": "#D92D20",
    "Komşu ARG (Direnç Geni)": "#F04438",
    "Transposon / IS Elementi": "#F79009",
    "İntegron / Rekombinaz": "#9E77ED",
    "Plazmid / Konjugasyon": "#2E90FA",
    "Bakteriyofaj (Prophage)": "#12B76A",
    "Hücresel / Diğer": "#98A2B3"
}

def plot_gene_neighborhood(neighbors_df: pd.DataFrame, target_arg: str, contig: str) -> go.Figure:
    fig = go.Figure()
    if neighbors_df.empty:
        fig.update_layout(title="Bu pencerede görüntülenecek gen bulunamadı.")
        return fig

    min_x = neighbors_df["start"].min()
    max_x = neighbors_df["end"].max()
    span = max(1, max_x - min_x)

    fig.add_shape(
        type="line",
        x0=min_x - span * 0.03,
        y0=0,
        x1=max_x + span * 0.03,
        y1=0,
        line=dict(color="#475467", width=3)
    )

    shown_legends = set()

    for _, row in neighbors_df.iterrows():
        s = row["start"]
        e = row["end"]
        strand = row["strand"]
        role = row["visual_role"]
        gene = row["gene_name"]
        product = row["product"]
        color = ROLE_COLORS.get(role, "#98A2B3")
        
        gene_len = e - s
        tip_len = min(gene_len * 0.35, span * 0.025)
        h_body = 0.22
        h_head = 0.36

        if strand == "+":
            neck = e - tip_len
            xs = [s, neck, neck, e, neck, neck, s, s]
            ys = [-h_body, -h_body, -h_head, 0, h_head, h_body, h_body, -h_body]
        else:
            neck = s + tip_len
            xs = [e, neck, neck, s, neck, neck, e, e]
            ys = [-h_body, -h_body, -h_head, 0, h_head, h_body, h_body, -h_body]

        show_leg = role not in shown_legends
        if show_leg:
            shown_legends.add(role)

        hover_txt = (
            f"<b>Gen:</b> {gene} ({strand})<br>"
            f"<b>Kategori:</b> {role}<br>"
            f"<b>Koordinat:</b> {s:,} – {e:,} bp ({gene_len:,} bp)<br>"
            f"<b>Ürün:</b> {product}"
        )

        fig.add_trace(go.Scatter(
            x=xs,
            y=ys,
            fill="toself",
            fillcolor=color,
            line=dict(color="#1D2939", width=1.2),
            mode="lines",
            name=role,
            legendgroup=role,
            showlegend=show_leg,
            text=hover_txt,
            hoverinfo="text"
        ))

        mid_x = (s + e) / 2
        label_y = 0.52 if strand == "+" else -0.55
        fig.add_annotation(
            x=mid_x,
            y=label_y,
            text=f"<b>{gene}</b>" if "ARG" in role else gene,
            showarrow=False,
            font=dict(size=11, color="#101828"),
            textangle=-25 if len(neighbors_df) > 8 else 0
        )

    fig.update_layout(
        title=dict(
            text=f"Genomik Bağlam Haritası: <b>{target_arg}</b> ({contig})",
            font=dict(size=16, color="#101828")
        ),
        xaxis=dict(title="Genomik Koordinat (bp)", showgrid=True, gridcolor="#EAECF0", zeroline=False),
        yaxis=dict(range=[-1.0, 1.0], visible=False),
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#FFFFFF",
        height=380,
        margin=dict(l=30, r=30, t=60, b=50),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    return fig

def plot_risk_bar_chart(summary_df: pd.DataFrame) -> go.Figure:
    df_sorted = summary_df.sort_values("HGT_Risk_Skoru", ascending=True)
    color_map = {
        "YÜKSEK RİSK (Kritik Mobilite)": "#D92D20",
        "ORTA RİSK (Potansiyel Mobil)": "#F79009",
        "DÜŞÜK RİSK (Kromozomal / Sabit)": "#12B76A"
    }
    fig = px.bar(
        df_sorted,
        x="HGT_Risk_Skoru",
        y="ARG_Geni",
        color="Risk_Seviyesi",
        color_discrete_map=color_map,
        orientation="h",
        text="HGT_Risk_Skoru",
        hover_data=["Contig", "Ilac_Sinifi", "Komsu_MGE_Sayisi"],
        title="Antibiyotik Direnç Genlerinin Yatay Transfer (HGT) Risk Sıralaması"
    )
    fig.update_layout(
        xaxis=dict(range=[0, 10.5], title="HGT Mobilite Risk Skoru (0 - 10)"),
        yaxis=dict(title="Direnç Geni (ARG)"),
        plot_bgcolor="#FFFFFF",
        height=max(300, len(df_sorted) * 48),
        legend_title="Risk Kategorisi"
    )
    fig.update_traces(textposition="outside")
    return fig
