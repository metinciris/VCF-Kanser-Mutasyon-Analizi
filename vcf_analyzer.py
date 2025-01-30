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
            # Kromozom numarasını düzenle
            chrom_num = chrom
            if chrom.upper() == 'X':
                chrom_num = '23'
            elif chrom.upper() == 'Y':
                chrom_num = '24'
            elif chrom.upper() == 'M':
                chrom_num = 'MT'
            
            # HGVS formatları
            hgvs_formats = [
                f"NC_0000{chrom_num if int(chrom_num) > 9 else '0'+chrom_num}.10:g.{pos}{ref}>{alt}" if chrom_num.isdigit() else f"NC_0000{chrom}.10:g.{pos}{ref}>{alt}",  # GRCh37
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
    def get_variant_effect(self, ref, alt):
        """Varyant tipini ve etkisini belirle"""
        if len(ref) == len(alt):
            if len(ref) == 1:
                return 'SNV'
            else:
                return 'MNV'
        elif len(ref) > len(alt):
            return 'Deletion'
        else:
            return 'Insertion'

    def analyze_variant(self, chrom, pos, ref, alt, variant_id=None):
        """Geliştirilmiş varyant analizi"""
        gene_info = self.get_gene_info_from_ucsc(chrom, pos)
        clinvar_info = self.get_clinvar_info(chrom, pos, ref, alt, variant_id)
        variant_type = self.get_variant_effect(ref, alt)

        # Varyant etkisini belirle
        effect = ''
        if len(ref) != len(alt):  # indel
            if len(ref) > len(alt):
                effect = 'frameshift_variant' if (len(ref) - len(alt)) % 3 != 0 else 'inframe_deletion'
            else:
                effect = 'frameshift_variant' if (len(alt) - len(ref)) % 3 != 0 else 'inframe_insertion'
        else:  # SNV
            if ref in ['A', 'G'] and alt in ['A', 'G']:  # transition
                effect = 'transition'
            elif ref in ['C', 'T'] and alt in ['C', 'T']:  # transition
                effect = 'transition'
            else:  # transversion
                effect = 'transversion'

        # Gene info kontrolü
        gene_symbol = gene_info.get('gene_symbol') if gene_info else ''
        is_cancer_gene = gene_info.get('is_cancer_gene', False) if gene_info else False
        transcript = gene_info.get('transcript', '') if gene_info else ''

        result = {
            'chromosome': chrom,
            'position': pos,
            'reference': ref,
            'alternate': alt,
            'variant_type': variant_type,
            'variant_effect': effect,
            'gene_symbol': gene_symbol,
            'is_cancer_gene': is_cancer_gene,
            'transcript': transcript,
            'in_clinvar': clinvar_info['found'],
            'clinvar_significance': clinvar_info['significance'],
            'clinvar_review_status': clinvar_info['review_status'],
            'clinvar_last_evaluated': clinvar_info['last_evaluated'],
            'clinvar_submission_count': clinvar_info['submission_count'],
            'cancer_gene_description': self.cancer_genes.get(gene_symbol, '')
        }
        return result

class VCFAnalyzer(VariantAnnotator):
    def __init__(self):
        super().__init__()
        self.root = tk.Tk()
        self.root.title("SDÜ NGS Panel Varyant Analizi")
        self.setup_gui()

    def setup_gui(self):
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Panel seçimi (opsiyonel)
        self.panel_var = tk.StringVar(value="all")
        panel_frame = ttk.LabelFrame(self.main_frame, text="Panel Tipi (Opsiyonel)", padding="5")
        panel_frame.pack(fill=tk.X, pady=5)
        
        ttk.Radiobutton(panel_frame, text="Tüm Genler", 
                       variable=self.panel_var, value="all").pack(side=tk.LEFT)
        ttk.Radiobutton(panel_frame, text="Solid Tümör Paneli", 
                       variable=self.panel_var, value="solid").pack(side=tk.LEFT)
        ttk.Radiobutton(panel_frame, text="Akciğer Paneli", 
                       variable=self.panel_var, value="lung").pack(side=tk.LEFT)
        
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
            status_text = f"Analiz başlıyor: {patient_id}"
            if panel_type != "all":
                status_text += f" ({panel_type.upper()} Panel)"
            self.update_status(status_text)
            
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
            
            # Panel raporu sadece panel seçiliyse oluşturulur
            if panel_type != "all":
                panel_report = self.create_panel_specific_report(df, panel_type)
                self.save_excel_report(df, patient_id, panel_type, panel_report)
            else:
                self.save_excel_report(df, patient_id, "all", None)
            
            self.create_visualizations(df, patient_id, panel_type)
            
            self.update_status(f"Analiz tamamlandı!", 100)
            
            completion_message = f"Toplam {len(df)} varyant analiz edildi."
            if panel_type != "all":
                panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                panel_variants = df[df['gene_symbol'].isin(panel_genes)]
                completion_message += f"\nPanel genlerinde {len(panel_variants)} varyant bulundu."
            messagebox.showinfo("Analiz Tamamlandı", completion_message)
            
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
        
        # 2. Gen dağılımı
        plt.subplot(2, 3, 2)
        if panel_type != "all":
            panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
            panel_df = df[df['gene_symbol'].isin(panel_genes)]
            if not panel_df.empty:
                sns.countplot(data=panel_df, y='gene_symbol', 
                            order=panel_df['gene_symbol'].value_counts().index)
                plt.title(f'{panel_type.upper()} Panel Genleri Dağılımı')
        else:
            cancer_df = df[df['is_cancer_gene']]
            if not cancer_df.empty:
                sns.countplot(data=cancer_df, y='gene_symbol',
                            order=cancer_df['gene_symbol'].value_counts().head(15).index)
                plt.title('En Sık Görülen Kanser Genleri (Top 15)')
        
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
        if not pathway_series.empty:
            sns.countplot(y=pathway_series, order=pathway_series.value_counts().head(10).index)
            plt.title('En Sık Görülen Yolaklar (Top 10)')
        
        # 5. Kanser geni oranı
        plt.subplot(2, 3, 5)
        if panel_type != "all":
            panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
            ratio_labels = ['Panel Genleri', 'Diğer Genler']
            ratio_values = [
                len(df[df['gene_symbol'].isin(panel_genes)]),
                len(df[~df['gene_symbol'].isin(panel_genes)])
            ]
        else:
            ratio_labels = ['Kanser Genleri', 'Diğer Genler']
            ratio_values = [
                len(df[df['is_cancer_gene']]),
                len(df[~df['is_cancer_gene']])
            ]
        plt.pie(ratio_values, labels=ratio_labels, autopct='%1.1f%%')
        plt.title('Gen Dağılım Oranı')
        
        plt.tight_layout()
        plt.savefig(f'{patient_id}_analysis_plots.png', dpi=300, bbox_inches='tight')
        plt.close()

    def save_excel_report(self, df, patient_id, panel_type, panel_report=None):
        """Geliştirilmiş Excel raporu"""
        try:
            output_file = f"{patient_id}_mutation_analysis.xlsx"
            with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
                workbook = writer.book
                
                # Tüm varyantlar
                df.to_excel(writer, sheet_name='All_Variants', index=False)
                
                # Kanser genleri
                cancer_variants = df[df['is_cancer_gene']].sort_values(
                    by=['variant_effect', 'gene_symbol'],
                    ascending=[False, True]
                )
                if not cancer_variants.empty:
                    cancer_variants.to_excel(writer, sheet_name='Cancer_Genes', index=False)
                
                # Panel-spesifik varyantlar
                if panel_type != "all":
                    panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                    panel_variants = df[df['gene_symbol'].isin(panel_genes)]
                    if not panel_variants.empty:
                        panel_variants.to_excel(writer, sheet_name='Panel_Variants', index=False)
                
                # Yüksek etkili varyantlar
                high_impact = df[df['variant_effect'].isin(['frameshift_variant', 'inframe_deletion'])]
                if not high_impact.empty:
                    high_impact.to_excel(writer, sheet_name='High_Impact_Variants', index=False)
                
                # Özet istatistikler
                summary_data = {
                    'Metrik': [
                        'Toplam Varyant Sayısı',
                        'Kanser Geni Varyantları',
                        'ClinVar Varyantları',
                        'Frameshift Varyant Sayısı',
                        'SNV Sayısı',
                        'İnsersiyon Sayısı',
                        'Delesyon Sayısı',
                        'Transition Sayısı',
                        'Transversion Sayısı'
                    ],
                    'Değer': [
                        len(df),
                        len(cancer_variants),
                        df['in_clinvar'].sum(),
                        len(df[df['variant_effect'] == 'frameshift_variant']),
                        len(df[df['variant_type'] == 'SNV']),
                        len(df[df['variant_type'] == 'Insertion']),
                        len(df[df['variant_type'] == 'Deletion']),
                        len(df[df['variant_effect'] == 'transition']),
                        len(df[df['variant_effect'] == 'transversion'])
                    ]
                }
                
                if panel_type != "all":
                    panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                    panel_variants = df[df['gene_symbol'].isin(panel_genes)]
                    summary_data['Metrik'].extend([
                        'Panel Genleri Varyantları',
                        'Panel Genleri Oranı (%)'
                    ])
                    summary_data['Değer'].extend([
                        len(panel_variants),
                        round(len(panel_variants) / len(df) * 100, 2)
                    ])
                
                pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
            
            print(f"Excel raporu başarıyla kaydedildi: {output_file}")
            # Excel dosyasını otomatik aç
            if os.path.exists(output_file):
                os.startfile(output_file)
                
        except Exception as e:
            print(f"Excel raporu oluşturulurken hata oluştu: {str(e)}")
            messagebox.showerror("Hata", f"Excel raporu oluşturulamadı: {str(e)}")

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
