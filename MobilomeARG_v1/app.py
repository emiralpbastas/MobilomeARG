import streamlit as st
import pandas as pd
import io
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
from parser import parse_gff3, parse_rgi, parse_plasmidfinder
from analyzer import analyze_hgt_risk
from visualizer import plot_gene_neighborhood, plot_risk_bar_chart

st.set_page_config(
    page_title="MobilomeARG | HGT & Genomik Bağlam Analizörü",
    page_icon="",
    layout="wide"
)

st.title("MobilomeARG: Genomik Bağlam & Yatay Gen Transferi (HGT) Risk Profilleyici")
st.markdown(
    "Short-read montajlarda (assembly) parçalanan plazmidler nedeniyle gözden kaçan **mobil antibiyotik direnç genlerini (ARG)**, "
    "genomik komşuluklarındaki **Transpozon (IS), İntegron ve Konjugasyon** elementlerini tarayarak tespit eder."
)

st.sidebar.header("Analiz Parametreleri")
window_kb = st.sidebar.slider(
    "Tarama Penceresi (± kb)",
    min_value=2,
    max_value=25,
    value=10,
    step=1,
    help="Her direnç geninin sağında ve solunda taranacak baz çifti mesafesi."
)
window_bp = window_kb * 1000

st.sidebar.markdown("---")
st.sidebar.subheader("Veri Yükleme")
rgi_file = st.sidebar.file_uploader("1. CARD / RGI Çıktısı (.txt / .tsv)", type=["txt", "tsv"])
gff_file = st.sidebar.file_uploader("2. Prokka / Bakta Anotasyonu (.gff3 / .gff)", type=["gff3", "gff"])
pf_file = st.sidebar.file_uploader("3. PlasmidFinder Çıktısı [Opsiyonel] (.tsv)", type=["tsv", "txt"])

if rgi_file is None or gff_file is None:
    st.info("Analizi başlatmak için sol menüden **CARD/RGI** (`.txt`/`.tsv`) ve **Prokka/Bakta** (`.gff3`/`.gff`) dosyalarınızı yükleyin.")
    st.stop()

try:
    gff_df = parse_gff3(gff_file)
    rgi_df = parse_rgi(rgi_file)
    plasmid_set = parse_plasmidfinder(pf_file)

    summary_df, neighborhoods = analyze_hgt_risk(gff_df, rgi_df, window_bp=window_bp, plasmid_contigs=plasmid_set)

    total_args = len(summary_df)
    high_risk = len(summary_df[summary_df["HGT_Risk_Skoru"] >= 7])
    mod_risk = len(summary_df[(summary_df["HGT_Risk_Skoru"] >= 3) & (summary_df["HGT_Risk_Skoru"] < 7)])
    low_risk = len(summary_df[summary_df["HGT_Risk_Skoru"] < 3])

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Toplam Direnç Geni (ARG)", total_args)
    c2.metric("Kritik Mobilite (7-10 Puan)", high_risk)
    c3.metric("Potansiyel Mobil (3-6 Puan)", mod_risk)
    c4.metric("Kromozomal / Sabit (0-2 Puan)", low_risk)

    st.markdown("---")
    st.subheader("İnteraktif Genomik Komşuluk Haritası (Gene Arrow Plot)")
    options = [f"{r['ARG_Geni']} ({r['Contig']}) — Skor: {r['HGT_Risk_Skoru']}/10" for _, r in summary_df.iterrows()]
    selected_idx = st.selectbox("İncelemek istediğiniz direnç genini seçin:", range(len(options)), format_func=lambda i: options[i])
    
    sel_row = summary_df.iloc[selected_idx]
    key = f"{sel_row['ARG_Geni']}_{sel_row['Contig']}"
    fig_arrow = plot_gene_neighborhood(neighborhoods[key], sel_row["ARG_Geni"], sel_row["Contig"])
    st.plotly_chart(fig_arrow, use_container_width=True)

    st.info(f"**Skor Gerekçesi:** {sel_row['Kanit_ve_Puan_Detayi']}")

    col_left, col_right = st.columns([1.1, 1.4])
    with col_left:
        st.subheader("HGT Risk Dağılımı")
        fig_bar = plot_risk_bar_chart(summary_df)
        st.plotly_chart(fig_bar, use_container_width=True)

    with col_right:
        st.subheader("Detaylı ARG & Mobilite Rapor Tablosu")
        st.dataframe(summary_df, use_container_width=True, height=360)

    st.markdown("---")
    st.subheader("Raporu Dışa Aktar")
    exp_col1, exp_col2 = st.columns(2)

    csv_bytes = summary_df.to_csv(index=False).encode("utf-8-sig")
    exp_col1.download_button(label="CSV Olarak İndir", data=csv_bytes, file_name="MobilomeARG_HGT_Raporu.csv", mime="text/csv")

    excel_buf = io.BytesIO()
    with pd.ExcelWriter(excel_buf, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="HGT_Risk_Ozeti", index=False)
    exp_col2.download_button(
        label="Excel (.xlsx) Olarak İndir",
        data=excel_buf.getvalue(),
        file_name="MobilomeARG_HGT_Raporu.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

except Exception as e:
    st.error(f"Dosya işlenirken bir hata oluştu: {e}")