import requests
import json
import time
from datetime import datetime
import pandas as pd
import tkinter as tk
from tkinter import filedialog, ttk
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

class VariantAnnotator:
    def __init__(self):
        self.cache = {}
        # NGS Tümör Paneli genleri
        self.cancer_genes = {
            'BRCA1': 'Herediter Meme/Over Kanseri',
            'BRCA2': 'Herediter Meme/Over Kanseri',
            'TP53': 'Li-Fraumeni Sendromu',
            'ALK': 'Akciğer Kanseri',
            'EGFR': 'Akciğer Kanseri',
            'KRAS': 'Kolorektal/Akciğer Kanseri',
            'BRAF': 'Melanom/Kolorektal Kanser',
            'NRAS': 'Melanom/Lösemi',
            'KIT': 'GIST',
            'PDGFRA': 'GIST',
            'IDH1': 'Glioma',
            'IDH2': 'Glioma',
            'FGFR1': 'Mesane Kanseri',
            'FGFR2': 'Mesane Kanseri',
            'FGFR3': 'Mesane Kanseri',
            'PIK3CA': 'Meme Kanseri',
            'AKT1': 'Meme Kanseri',
            'ERBB2': 'Meme Kanseri (HER2)',
            'PTEN': 'Cowden Sendromu',
            'CDH1': 'Herediter Diffüz Mide Kanseri',
            'BMPR1A': 'Juvenil Polipozis',
            'SMAD4': 'Juvenil Polipozis',
            'STK11': 'Peutz-Jeghers Sendromu',
            'MLH1': 'Lynch Sendromu',
            'MSH2': 'Lynch Sendromu',
            'MSH6': 'Lynch Sendromu',
            'PMS2': 'Lynch Sendromu',
            'EPCAM': 'Lynch Sendromu',
            'MPL': 'Trombositemi',
            'TPM3': 'Papiller Tiroid Kanseri',
            'NTRK1': 'Tiroid ve Diğer Kanserler',
            'SLC45A3': 'Prostat Kanseri Füzyon Partner',
            'ARID1A': 'Over ve Endometrium Kanseri'
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
                    'gene_type': None
                }

                if 'refGene' in data and data['refGene']:
                    ref_gene = data['refGene'][0]
                    gene_info.update({
                        'gene_symbol': ref_gene.get('name2'),
                        'transcript': ref_gene.get('name'),
                        'strand': ref_gene.get('strand'),
                        'exon_count': ref_gene.get('exonCount')
                    })
                    gene_info['is_cancer_gene'] = gene_info['gene_symbol'] in self.cancer_genes

                self.cache[cache_key] = gene_info
                return gene_info
        except Exception as e:
            print(f"UCSC API hatası: {str(e)}")
        return None

    def get_clinvar_info(self, chrom, pos, ref, alt, variant_id=None):
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
                
                # Klinik önem bilgisini al
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

    def get_impact_level(self, mutation_type):
        """Mutasyon etkisini belirle"""
        for impact, types in self.impact_levels.items():
            if any(t in mutation_type.lower() for t in types):
                return impact
        return 'UNKNOWN'

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

        result = {
            'chromosome': chrom,
            'position': pos,
            'reference': ref,
            'alternate': alt,
            'variant_type': variant_type,
            'variant_effect': effect,
            'gene_symbol': gene_info['gene_symbol'] if gene_info else '',
            'is_cancer_gene': gene_info['is_cancer_gene'] if gene_info else False,
            'transcript': gene_info['transcript'] if gene_info else '',
            'in_clinvar': clinvar_info['found'],
            'clinvar_significance': clinvar_info['significance'],
            'clinvar_review_status': clinvar_info['review_status'],
            'clinvar_last_evaluated': clinvar_info['last_evaluated'],
            'clinvar_submission_count': clinvar_info['submission_count'],
            'cancer_gene_description': self.cancer_genes.get(gene_info['gene_symbol'] if gene_info else None, '')
        }
        return result

    def create_summary_report(self, df, cancer_variants):
        """Detaylı özet raporu oluştur"""
        high_impact_variants = df[df['variant_effect'].isin(['frameshift_variant', 'inframe_deletion'])]
        
        summary = {
            'Genel Özet': {
                'Toplam Varyant': len(df),
                'Kanser Geni Varyantları': len(cancer_variants),
                'Yüksek Etkili Varyantlar': len(high_impact_variants)
            },
            'Önemli Bulgular': [],
            'Kanser Genleri Analizi': []
        }
        
        # Önemli bulguları belirle
        for _, variant in high_impact_variants.iterrows():
            if variant['is_cancer_gene']:
                summary['Önemli Bulgular'].append(
                    f"{variant['gene_symbol']}: {variant['variant_effect']} "
                    f"(chr{variant['chromosome']}:{variant['position']} "
                    f"{variant['reference']}>{variant['alternate']})"
                )
        
        # Her kanser geni için detaylı analiz
        for gene in cancer_variants['gene_symbol'].unique():
            gene_variants = cancer_variants[cancer_variants['gene_symbol'] == gene]
            summary['Kanser Genleri Analizi'].append({
                'Gen': gene,
                'Kanser Tipi': self.cancer_genes[gene],
                'Varyant Sayısı': len(gene_variants),
                'Varyant Tipleri': ', '.join(gene_variants['variant_effect'].unique())
            })
        
        return summary
      class VCFAnalyzer(VariantAnnotator):
    def __init__(self):
        super().__init__()
        self.root = tk.Tk()
        self.root.title("VCF Kanser Mutasyon Analizi")
        self.setup_gui()

    def setup_gui(self):
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
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
            self.update_status(f"Analiz başlıyor: {patient_id}")
            
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
            cancer_variants = df[df['is_cancer_gene']]
            
            summary_report = self.create_summary_report(df, cancer_variants)
            self.save_excel_report(df, patient_id, summary_report)
            self.create_visualizations(df, patient_id)
            
            self.update_status(f"Analiz tamamlandı!", 100)
            
        except Exception as e:
            self.update_status(f"Hata: {str(e)}")

    def create_visualizations(self, df, patient_id):
        """Geliştirilmiş görselleştirmeler"""
        sns.set_theme(style="whitegrid")
        
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Varyant tipi dağılımı
        plt.subplot(2, 3, 1)
        sns.countplot(data=df, x='variant_type')
        plt.title('Varyant Tipi Dağılımı')
        plt.xticks(rotation=45)
        
        # 2. Kanser genleri
        plt.subplot(2, 3, 2)
        cancer_df = df[df['is_cancer_gene']]
        if not cancer_df.empty:
            sns.countplot(data=cancer_df, y='gene_symbol', 
                         order=cancer_df['gene_symbol'].value_counts().index)
        plt.title('Kanser Genleri Dağılımı')
        
        # 3. Kromozom dağılımı
        plt.subplot(2, 3, 3)
        sns.countplot(data=df, x='chromosome')
        plt.title('Kromozom Dağılımı')
        
        # 4. Varyant etki dağılımı
        plt.subplot(2, 3, 4)
        sns.countplot(data=df, y='variant_effect')
        plt.title('Varyant Etki Dağılımı')
        
        # 5. Kanser geni vs normal gen oranı
        plt.subplot(2, 3, 5)
        cancer_ratio = pd.Series({
            'Kanser Genleri': len(df[df['is_cancer_gene']]),
            'Diğer Genler': len(df[~df['is_cancer_gene']])
        })
        plt.pie(cancer_ratio, labels=cancer_ratio.index, autopct='%1.1f%%')
        plt.title('Kanser Geni Oranı')
        
        plt.tight_layout()
        plt.savefig(f'{patient_id}_analysis_plots.png', dpi=300, bbox_inches='tight')
        plt.close()

    def save_excel_report(self, df, patient_id, summary_report):
        """Geliştirilmiş Excel raporu"""
        output_file = f"{patient_id}_mutation_analysis.xlsx"
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            workbook = writer.book
            
            # Formatlar
            header_format = workbook.add_format({
                'bold': True,
                'bg_color': '#D7E4BC',
                'border': 1
            })
            warning_format = workbook.add_format({
                'bg_color': '#FFC7CE',
                'font_color': '#9C0006'
            })
            
            # Tüm varyantlar
            df.to_excel(writer, sheet_name='All_Variants', index=False)
            
            # Kanser genleri - önem sırasına göre sırala
            cancer_variants = df[df['is_cancer_gene']].sort_values(
                by=['variant_effect', 'gene_symbol'],
                ascending=[False, True]
            )
            cancer_variants.to_excel(writer, sheet_name='Cancer_Genes', index=False)
            
            # Yüksek etkili varyantlar
            high_impact = df[
                (df['variant_effect'].isin(['frameshift_variant', 'inframe_deletion'])) &
                (df['is_cancer_gene'])
            ]
            high_impact.to_excel(writer, sheet_name='High_Impact_Variants', index=False)
            
            # Özet raporu
            summary_df = pd.DataFrame([
                {'Metrik': k, 'Değer': v} for k, v in summary_report['Genel Özet'].items()
            ])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Kanser genleri özeti
            cancer_summary = pd.DataFrame(summary_report['Kanser Genleri Analizi'])
            cancer_summary.to_excel(writer, sheet_name='Cancer_Genes_Summary', index=False)
            
            # Önemli bulgular
            if summary_report['Önemli Bulgular']:
                findings_df = pd.DataFrame({'Önemli Bulgular': summary_report['Önemli Bulgular']})
                findings_df.to_excel(writer, sheet_name='Important_Findings', index=False)

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
