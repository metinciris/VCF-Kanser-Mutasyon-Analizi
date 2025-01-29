import requests
import json
import time
from datetime import datetime
import pandas as pd
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

class VariantAnnotator:
    def __init__(self):
        self.cache = {}
        # SDÜ Solid Tümör ve Lung Paneli genleri
        self.cancer_genes = {
            'AKT1': 'PI3K/AKT Yolağı',
            'AKT2': 'PI3K/AKT Yolağı',
            'AKT3': 'PI3K/AKT Yolağı',
            'ALK': 'Akciğer Kanseri',
            'APC': 'Kolorektal Kanser',
            'ARID1A': 'Kromatin Düzenleyici',
            'ATM': 'DNA Hasar Yanıtı',
            'ATRX': 'Kromatin Düzenleyici',
            'AURKA': 'Hücre Döngüsü',
            'BAP1': 'Tümör Supresör',
            'BRAF': 'MAPK Yolağı',
            'BRCA1': 'DNA Onarımı',
            'BRCA2': 'DNA Onarımı',
            'CCND1': 'Hücre Döngüsü',
            'CCNE1': 'Hücre Döngüsü',
            'CDK12': 'Hücre Döngüsü',
            'CDK4': 'Hücre Döngüsü',
            'CDKN2A': 'Hücre Döngüsü',
            'CHEK1': 'DNA Hasar Yanıtı',
            'CHEK2': 'DNA Hasar Yanıtı',
            'CREBBP': 'Transkripsiyon Düzenleyici',
            'CSF1R': 'Tirozin Kinaz',
            'CTNNB1': 'WNT Yolağı',
            'DDR2': 'Tirozin Kinaz',
            'DNMT3A': 'Epigenetik Düzenleyici',
            'EGFR': 'Tirozin Kinaz',
            'EPCAM': 'Hücre Adezyonu',
            'ERBB2': 'HER2',
            'ERBB3': 'HER Ailesi',
            'ERBB4': 'HER Ailesi',
            'ESR1': 'Hormon Reseptörü',
            'FAT1': 'WNT Yolağı',
            'FBXO11': 'Protein Degradasyonu',
            'FGFR1': 'Tirozin Kinaz',
            'FGFR2': 'Tirozin Kinaz',
            'FGFR3': 'Tirozin Kinaz',
            'FGFR4': 'Tirozin Kinaz',
            'FLT3': 'Tirozin Kinaz',
            'GNA11': 'G Protein',
            'GNAQ': 'G Protein',
            'GNAS': 'G Protein',
            'HRAS': 'RAS Yolağı',
            'IDH1': 'Metabolizma',
            'IDH2': 'Metabolizma',
            'KEAP1': 'Oksidatif Stres',
            'KIT': 'Tirozin Kinaz',
            'KRAS': 'RAS Yolağı',
            'MAP2K1': 'MAPK Yolağı',
            'MAP2K2': 'MAPK Yolağı',
            'MDM2': 'p53 Yolağı',
            'MET': 'Tirozin Kinaz',
            'MLH1': 'DNA Onarımı',
            'MPL': 'Trombopoetin Reseptörü',
            'MSH2': 'DNA Onarımı',
            'MSH6': 'DNA Onarımı',
            'MTOR': 'PI3K/AKT/mTOR Yolağı',
            'MYC': 'Transkripsiyon Faktörü',
            'NF1': 'RAS Yolağı',
            'NOTCH1': 'NOTCH Yolağı',
            'NOTCH2': 'NOTCH Yolağı',
            'NOTCH3': 'NOTCH Yolağı',
            'NOTCH4': 'NOTCH Yolağı',
            'NRAS': 'RAS Yolağı',
            'NTRK1': 'Tirozin Kinaz',
            'NTRK2': 'Tirozin Kinaz',
            'NTRK3': 'Tirozin Kinaz',
            'PDGFRA': 'Tirozin Kinaz',
            'PIK3CA': 'PI3K Yolağı',
            'PIK3R1': 'PI3K Yolağı',
            'POLE': 'DNA Polimeraz',
            'PTEN': 'PI3K Yolağı',
            'PTPN11': 'Tirozin Fosfataz',
            'RAF1': 'MAPK Yolağı',
            'RB1': 'Hücre Döngüsü',
            'RET': 'Tirozin Kinaz',
            'ROS1': 'Tirozin Kinaz',
            'SETD2': 'Histon Metiltransferaz',
            'SMAD4': 'TGF-beta Yolağı',
            'SRC': 'Tirozin Kinaz',
            'STK11': 'Tümör Supresör',
            'TERT': 'Telomeraz',
            'TET2': 'Epigenetik Düzenleyici',
            'TP53': 'Tümör Supresör',
            'TSC1': 'mTOR Yolağı',
            'TSC2': 'mTOR Yolağı',
            'VHL': 'Hipoksi Yolağı'
        }

        # Panel-spesifik genler
        self.solid_panel_genes = {
            'BRCA1', 'BRCA2', 'TP53', 'APC', 'KRAS', 'PIK3CA', 'PTEN', 'MLH1', 
            'MSH2', 'MSH6', 'EPCAM', 'ATM', 'CHEK2', 'BRAF', 'EGFR'
        }
        
        self.lung_panel_genes = {
            'EGFR', 'ALK', 'ROS1', 'KRAS', 'BRAF', 'MET', 'ERBB2', 'RET', 
            'NTRK1', 'NTRK2', 'NTRK3', 'PIK3CA', 'STK11'
        }

        # Varyant etki seviyeleri
        self.impact_levels = {
            'HIGH': ['frameshift', 'nonsense', 'splice_acceptor', 'splice_donor', 'stop_gained', 'stop_lost'],
            'MODERATE': ['missense', 'inframe_insertion', 'inframe_deletion', 'protein_altering'],
            'LOW': ['synonymous', 'splice_region', 'intronic', 'upstream', 'downstream'],
            'MODIFIER': ['intergenic', 'non_coding_transcript', 'regulatory_region']
        }
    def get_gene_info_from_ucsc(self, chrom, pos):
        """UCSC API'den geliştirilmiş gen bilgisi"""
        cache_key = f"ucsc_{chrom}:{pos}"
        if cache_key in self.cache:
            return self.cache[cache_key]

        url = f"https://api.genome.ucsc.edu/getData/track?"
        params = {
            'genome': 'hg19',
            'track': 'refGene,knownGene,ensGene',
            'chrom': f'chr{chrom}',
            'start': int(pos)-1,
            'end': pos
        }

        try:
            response = requests.get(url, params=params, timeout=10)
            if response.ok:
                data = response.json()
                gene_info = {
                    'gene_symbol': None,
                    'transcript': None,
                    'strand': None,
                    'exon_count': None,
                    'is_cancer_gene': False,
                    'gene_type': None,
                    'in_solid_panel': False,
                    'in_lung_panel': False
                }

                if 'refGene' in data and data['refGene']:
                    ref_gene = data['refGene'][0]
                    gene_symbol = ref_gene.get('name2')
                    gene_info.update({
                        'gene_symbol': gene_symbol,
                        'transcript': ref_gene.get('name'),
                        'strand': ref_gene.get('strand'),
                        'exon_count': ref_gene.get('exonCount'),
                        'is_cancer_gene': gene_symbol in self.cancer_genes,
                        'in_solid_panel': gene_symbol in self.solid_panel_genes,
                        'in_lung_panel': gene_symbol in self.lung_panel_genes
                    })

                self.cache[cache_key] = gene_info
                return gene_info
        except Exception as e:
            print(f"UCSC API hatası: {str(e)}")
        return None

    def get_clinvar_info(self, chrom, pos, ref, alt, variant_id=None):
        """Geliştirilmiş ClinVar bilgisi"""
        try:
            # HGVS formatları
            hgvs_formats = [
                f"NC_0000{chrom if int(chrom) > 9 else '0'+chrom}.10:g.{pos}{ref}>{alt}",  # GRCh37
                f"chr{chrom}:g.{pos}{ref}>{alt}",  # Genomik
                f"{chrom}:{pos}{ref}>{alt}"  # Basit format
            ]
            
            for hgvs in hgvs_formats:
                params = {
                    'db': 'clinvar',
                    'term': f'"{hgvs}"[Variant Name]',
                    'retmode': 'json'
                }
                response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params=params)
                if response.ok:
                    data = response.json()
                    if int(data['esearchresult'].get('count', 0)) > 0:
                        return self._get_variant_details(data['esearchresult']['idlist'][0])

            # Genişletilmiş koordinat sorgusu
            params = {
                'db': 'clinvar',
                'term': f"{chrom}[Chr] AND {pos}[Base Position] AND {ref}[Reference allele] AND {alt}[Alternate allele]",
                'retmode': 'json'
            }
            response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params=params)
            if response.ok:
                data = response.json()
                if int(data['esearchresult'].get('count', 0)) > 0:
                    return self._get_variant_details(data['esearchresult']['idlist'][0])

        except Exception as e:
            print(f"ClinVar sorgu hatası: {str(e)}")
        
        return {
            'found': False,
            'clinvar_ids': [],
            'significance': '',
            'review_status': '',
            'last_evaluated': '',
            'submission_count': 0
        }

    def _get_variant_details(self, clinvar_id):
        """Geliştirilmiş ClinVar varyant detayları"""
        try:
            summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            summary_params = {
                'db': 'clinvar',
                'id': clinvar_id,
                'retmode': 'json'
            }
            
            response = requests.get(summary_url, params=summary_params)
            if response.ok:
                data = response.json()
                variant_info = data['result'][clinvar_id]
                
                clinical_significance = variant_info.get('clinical_significance', {})
                if isinstance(clinical_significance, dict):
                    significance = clinical_significance.get('description', '')
                else:
                    significance = str(clinical_significance)

                return {
                    'found': True,
                    'clinvar_ids': [clinvar_id],
                    'significance': significance,
                    'review_status': variant_info.get('review_status', ''),
                    'last_evaluated': variant_info.get('last_evaluated', ''),
                    'submission_count': variant_info.get('submission_count', 0)
                }
        except Exception as e:
            print(f"Varyant detay hatası: {str(e)}")
        
        return {
            'found': False,
            'clinvar_ids': [],
            'significance': '',
            'review_status': '',
            'last_evaluated': '',
            'submission_count': 0
        }

    def create_panel_specific_report(self, df, panel_type="solid"):
        """Panel tipine özel rapor oluştur"""
        panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
        panel_variants = df[df['gene_symbol'].isin(panel_genes)]
        
        report = {
            'Panel Tipi': panel_type.upper(),
            'Analiz Edilen Panel Genleri': len(panel_genes),
            'Bulunan Panel Varyantları': len(panel_variants),
            'Gen Bazlı Bulgular': {},
            'Yolak Analizi': {}
        }
        
        # Gen bazlı analiz
        for gene in panel_genes:
            gene_variants = panel_variants[panel_variants['gene_symbol'] == gene]
            if not gene_variants.empty:
                report['Gen Bazlı Bulgular'][gene] = {
                    'Varyant Sayısı': len(gene_variants),
                    'Varyant Tipleri': gene_variants['variant_effect'].unique().tolist(),
                    'Yolak': self.cancer_genes.get(gene, 'Bilinmiyor')
                }
        
        # Yolak analizi
        pathway_variants = {}
        for _, variant in panel_variants.iterrows():
            gene = variant['gene_symbol']
            pathway = self.cancer_genes.get(gene, 'Diğer')
            if pathway not in pathway_variants:
                pathway_variants[pathway] = 0
            pathway_variants[pathway] += 1
        
        report['Yolak Analizi'] = pathway_variants
        
        return report
class VCFAnalyzer(VariantAnnotator):
    def __init__(self):
        super().__init__()
        self.root = tk.Tk()
        self.root.title("SDÜ NGS Panel Varyant Analizi")
        self.setup_gui()

    def setup_gui(self):
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Panel seçimi
        self.panel_var = tk.StringVar(value="solid")
        ttk.Label(self.main_frame, text="Panel Tipi:").pack(pady=5)
        ttk.Radiobutton(self.main_frame, text="Solid Tümör Paneli", 
                       variable=self.panel_var, value="solid").pack()
        ttk.Radiobutton(self.main_frame, text="Akciğer Paneli", 
                       variable=self.panel_var, value="lung").pack()
        
        ttk.Button(self.main_frame, text="VCF Dosyası Seç",
                  command=self.select_file).pack(pady=5)
        
        self.progress = ttk.Progressbar(self.main_frame, length=300, mode='determinate')
        self.progress.pack(pady=5)
        
        self.status_var = tk.StringVar(value="Hazır")
        ttk.Label(self.main_frame, textvariable=self.status_var).pack(pady=5)

    def select_file(self):
        vcf_file = filedialog.askopenfilename(
            title="VCF Dosyası Seç",
            filetypes=[("VCF files", "*.vcf")]
        )
        if vcf_file:
            self.analyze_vcf(vcf_file)

    def update_status(self, message, progress=None):
        self.status_var.set(message)
        if progress is not None:
            self.progress['value'] = progress
        self.root.update()

    def analyze_vcf(self, vcf_file):
        try:
            patient_id = self.extract_patient_id(vcf_file)
            panel_type = self.panel_var.get()
            self.update_status(f"Analiz başlıyor: {patient_id} ({panel_type.upper()} Panel)")
            
            variants = self.read_vcf(vcf_file)
            total_variants = len(variants)
            results = []

            for i, variant in enumerate(variants):
                try:
                    result = self.analyze_variant(
                        str(variant['CHROM']),
                        str(variant['POS']),
                        str(variant['REF']),
                        str(variant['ALT']),
                        variant['ID']
                    )
                    results.append(result)
                    
                    if i % 10 == 0:
                        progress = (i + 1) * 100 / total_variants
                        self.update_status(
                            f"İşlenen: {i+1}/{total_variants} varyant",
                            progress
                        )
                except Exception as e:
                    print(f"Varyant analiz hatası: {str(e)}")

            df = pd.DataFrame(results)
            panel_report = self.create_panel_specific_report(df, panel_type)
            self.save_excel_report(df, patient_id, panel_report)
            self.create_visualizations(df, patient_id, panel_type)
            
            self.update_status(f"Analiz tamamlandı!", 100)
            messagebox.showinfo("Analiz Tamamlandı", 
                              f"Toplam {len(df)} varyant analiz edildi.\n"
                              f"Panel genlerinde {panel_report['Bulunan Panel Varyantları']} varyant bulundu.")
            
        except Exception as e:
            self.update_status(f"Hata: {str(e)}")
            messagebox.showerror("Hata", str(e))

    def create_visualizations(self, df, patient_id, panel_type):
        """Geliştirilmiş görselleştirmeler"""
        sns.set_theme(style="whitegrid")
        
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Varyant tipi dağılımı
        plt.subplot(2, 3, 1)
        sns.countplot(data=df, x='variant_type')
        plt.title('Varyant Tipi Dağılımı')
        plt.xticks(rotation=45)
        
        # 2. Panel genleri
        plt.subplot(2, 3, 2)
        panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
        panel_df = df[df['gene_symbol'].isin(panel_genes)]
        if not panel_df.empty:
            sns.countplot(data=panel_df, y='gene_symbol', 
                         order=panel_df['gene_symbol'].value_counts().index)
        plt.title(f'{panel_type.upper()} Panel Genleri Dağılımı')
        
        # 3. Varyant etki dağılımı
        plt.subplot(2, 3, 3)
        sns.countplot(data=df, y='variant_effect')
        plt.title('Varyant Etki Dağılımı')
        
        # 4. Yolak dağılımı
        plt.subplot(2, 3, 4)
        pathway_data = []
        for _, row in df[df['is_cancer_gene']].iterrows():
            pathway_data.append(self.cancer_genes.get(row['gene_symbol'], 'Diğer'))
        pathway_series = pd.Series(pathway_data)
        sns.countplot(y=pathway_series)
        plt.title('Yolak Dağılımı')
        
        # 5. Panel vs Non-panel gen oranı
        plt.subplot(2, 3, 5)
        panel_ratio = pd.Series({
            'Panel Genleri': len(df[df['gene_symbol'].isin(panel_genes)]),
            'Diğer Genler': len(df[~df['gene_symbol'].isin(panel_genes)])
        })
        plt.pie(panel_ratio, labels=panel_ratio.index, autopct='%1.1f%%')
        plt.title('Panel Genleri Oranı')
        
        plt.tight_layout()
        plt.savefig(f'{patient_id}_{panel_type}_panel_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

    def save_excel_report(self, df, patient_id, panel_report):
        """Geliştirilmiş Excel raporu"""
        output_file = f"{patient_id}_panel_analysis.xlsx"
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            workbook = writer.book
            
            # Formatlar
            header_format = workbook.add_format({
                'bold': True,
                'bg_color': '#D7E4BC',
                'border': 1
            })
            
            # Panel özet sayfası
            panel_summary = pd.DataFrame([
                {'Metrik': k, 'Değer': v} 
                for k, v in panel_report.items() 
                if not isinstance(v, dict)
            ])
            panel_summary.to_excel(writer, sheet_name='Panel_Summary', index=False)
            
            # Gen bazlı bulgular
            gene_findings = pd.DataFrame([
                {
                    'Gen': gene,
                    'Varyant Sayısı': info['Varyant Sayısı'],
                    'Varyant Tipleri': ', '.join(info['Varyant Tipleri']),
                    'Yolak': info['Yolak']
                }
                for gene, info in panel_report['Gen Bazlı Bulgular'].items()
            ])
            if not gene_findings.empty:
                gene_findings.to_excel(writer, sheet_name='Gene_Findings', index=False)
            
            # Yolak analizi
            pathway_analysis = pd.DataFrame([
                {'Yolak': k, 'Varyant Sayısı': v}
                for k, v in panel_report['Yolak Analizi'].items()
            ])
            if not pathway_analysis.empty:
                pathway_analysis.to_excel(writer, sheet_name='Pathway_Analysis', index=False)
            
            # Tüm varyantlar
            df.to_excel(writer, sheet_name='All_Variants', index=False)
            
            # Panel varyantları
            panel_genes = self.solid_panel_genes if panel_report['Panel Tipi'] == 'SOLID' else self.lung_panel_genes
            panel_variants = df[df['gene_symbol'].isin(panel_genes)]
            if not panel_variants.empty:
                panel_variants.to_excel(writer, sheet_name='Panel_Variants', index=False)

    def read_vcf(self, vcf_file):
        variants = []
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    variants.append({
                        'CHROM': fields[0],
                        'POS': fields[1],
                        'ID': fields[2],
                        'REF': fields[3],
                        'ALT': fields[4]
                    })
        return variants

    def extract_patient_id(self, vcf_path):
        try:
            with open(vcf_path, 'r') as f:
                for line in f:
                    if line.startswith('##fileOrigin='):
                        match = re.search(r'MP[0-9-]+', line)
                        if match:
                            return match.group(0)
            return os.path.basename(vcf_path).split('_')[0]
        except:
            return "Unknown"

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    analyzer = VCFAnalyzer()
    analyzer.run()
