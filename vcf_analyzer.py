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

class MultiDBVariantAnnotator:
    def __init__(self):
        self.cache = {}
        # API kaynakları
        self.api_sources = {
            "ClinVar": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            "gnomAD": "https://gnomad.broadinstitute.org/api/variant/",
            "MyVariant": "https://myvariant.info/v1/variant/",
            "dbSNP": "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
        }
        
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
                    'Gen': None,
                    'Transkript': None,
                    'Strand': None,
                    'Ekzon Sayısı': None,
                    'Kanser Geni': False,
                    'Gen Tipi': None,
                    'Solid Panel': False,
                    'Akciğer Panel': False
                }

                if 'refGene' in data and data['refGene']:
                    ref_gene = data['refGene'][0]
                    gene_symbol = ref_gene.get('name2')
                    gene_info.update({
                        'Gen': gene_symbol,
                        'Transkript': ref_gene.get('name'),
                        'Strand': ref_gene.get('strand'),
                        'Ekzon Sayısı': ref_gene.get('exonCount'),
                        'Kanser Geni': gene_symbol in self.cancer_genes,
                        'Solid Panel': gene_symbol in self.solid_panel_genes,
                        'Akciğer Panel': gene_symbol in self.lung_panel_genes
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
                f"NC_0000{chrom_num if int(chrom_num) > 9 else '0'+chrom_num}.10:g.{pos}{ref}>{alt}" if chrom_num.isdigit() else f"NC_0000{chrom}.10:g.{pos}{ref}>{alt}",
                f"chr{chrom}:g.{pos}{ref}>{alt}",
                f"{chrom}:{pos}{ref}>{alt}"
            ]
            
            for hgvs in hgvs_formats:
                params = {
                    'db': 'clinvar',
                    'term': f'"{hgvs}"[Variant Name]',
                    'retmode': 'json'
                }
                response = requests.get(self.api_sources['ClinVar'], params=params)
                if response.ok:
                    data = response.json()
                    if int(data['esearchresult'].get('count', 0)) > 0:
                        return self._get_variant_details(data['esearchresult']['idlist'][0])

        except Exception as e:
            print(f"ClinVar sorgu hatası: {str(e)}")
        
        return {
            'Bulundu': False,
            'ClinVar ID': [],
            'Klinik Önemi': '',
            'İnceleme Durumu': '',
            'Son Değerlendirme': '',
            'Başvuru Sayısı': 0
        }

    def get_gnomad_info(self, chrom, pos, ref, alt):
        """gnomAD bilgisi"""
        try:
            response = requests.get(f"{self.api_sources['gnomAD']}{chrom}-{pos}-{ref}-{alt}")
            if response.ok:
                data = response.json()
                return {
                    'Allel Frekansı': data.get('af', {}).get('af', ''),
                    'Filtre': data.get('filters', [])
                }
        except Exception as e:
            print(f"gnomAD sorgu hatası: {str(e)}")
        return {}

    def get_myvariant_info(self, chrom, pos, ref, alt):
        """MyVariant.info bilgisi"""
        try:
            hgvs = f"chr{chrom}:g.{pos}{ref}>{alt}"
            response = requests.get(f"{self.api_sources['MyVariant']}{hgvs}")
            if response.ok:
                data = response.json()
                return {
                    'CADD Skoru': data.get('cadd', {}).get('phred', ''),
                    'SIFT Tahmini': data.get('dbnsfp', {}).get('sift_pred', ''),
                    'PolyPhen Tahmini': data.get('dbnsfp', {}).get('polyphen2_hdiv_pred', '')
                }
        except Exception as e:
            print(f"MyVariant.info sorgu hatası: {str(e)}")
        return {}

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
                    'Bulundu': True,
                    'ClinVar ID': [clinvar_id],
                    'Klinik Önemi': significance,
                    'İnceleme Durumu': variant_info.get('review_status', ''),
                    'Son Değerlendirme': variant_info.get('last_evaluated', ''),
                    'Başvuru Sayısı': variant_info.get('submission_count', 0)
                }
        except Exception as e:
            print(f"Varyant detay hatası: {str(e)}")
        
        return {
            'Bulundu': False,
            'ClinVar ID': [],
            'Klinik Önemi': '',
            'İnceleme Durumu': '',
            'Son Değerlendirme': '',
            'Başvuru Sayısı': 0
        }

    def analyze_variant(self, chrom, pos, ref, alt, variant_id=None, dp='-', vaf='-'):
        """Geliştirilmiş varyant analizi"""
        gene_info = self.get_gene_info_from_ucsc(chrom, pos)
        
        # Veritabanı sorgularını kontrol et
        variant_annotations = {
            'ClinVar': {'Bulundu': False},
            'gnomAD': {},
            'MyVariant': {}
        }
        
        if hasattr(self, 'clinvar_var') and self.clinvar_var.get():
            clinvar_info = self.get_clinvar_info(chrom, pos, ref, alt, variant_id)
            variant_annotations['ClinVar'] = clinvar_info
        
        if hasattr(self, 'gnomad_var') and self.gnomad_var.get():
            gnomad_info = self.get_gnomad_info(chrom, pos, ref, alt)
            variant_annotations['gnomAD'] = gnomad_info
        
        if hasattr(self, 'myvariant_var') and self.myvariant_var.get():
            myvariant_info = self.get_myvariant_info(chrom, pos, ref, alt)
            variant_annotations['MyVariant'] = myvariant_info

        # Varyant tipini belirle
        if len(ref) == len(alt):
            variant_type = 'SNV' if len(ref) == 1 else 'MNV'
        elif len(ref) > len(alt):
            variant_type = 'Delesyon'
        else:
            variant_type = 'İnsersiyon'

        # Varyant etkisini belirle
        if len(ref) != len(alt):  # indel
            if len(ref) > len(alt):
                effect = 'frameshift_variant' if (len(ref) - len(alt)) % 3 != 0 else 'inframe_deletion'
                protein_impact = 'Protein Yapısını Bozan' if (len(ref) - len(alt)) % 3 != 0 else 'Protein Modifikasyonu'
            else:
                effect = 'frameshift_variant' if (len(alt) - len(ref)) % 3 != 0 else 'inframe_insertion'
                protein_impact = 'Protein Yapısını Bozan' if (len(alt) - len(ref)) % 3 != 0 else 'Protein Modifikasyonu'
            mutation_type = 'Fonksiyon Kaybı' if 'frameshift' in effect else 'Orta Etki'
        else:  # SNV
            if ref in ['A', 'G'] and alt in ['A', 'G']:
                effect = 'transition'
                protein_impact = 'A>G Değişimi'
            elif ref in ['C', 'T'] and alt in ['C', 'T']:
                effect = 'transition'
                protein_impact = 'C>T Değişimi'
            else:
                effect = 'transversion'
                protein_impact = f'{ref}>{alt} Değişimi'
            mutation_type = 'Değişim'

        # Protein fonksiyon tahmini
        if effect == 'frameshift_variant':
            protein_function = 'Muhtemel Fonksiyon Kaybı'
        elif effect in ['inframe_deletion', 'inframe_insertion']:
            protein_function = 'Kısmi Fonksiyon Değişimi'
        elif effect in ['transition', 'transversion']:
            protein_function = 'Olası Fonksiyon Değişimi' if gene_info and gene_info.get('Kanser Geni') else 'Belirsiz'
        else:
            protein_function = 'Belirsiz'

        # Klinik etki seviyesi
        if effect == 'frameshift_variant':
            clinical_impact = 'Yüksek'
        elif effect in ['inframe_deletion', 'inframe_insertion']:
            clinical_impact = 'Orta'
        elif effect in ['transition', 'transversion']:
            clinical_impact = 'Düşük'
        else:
            clinical_impact = 'Belirsiz'

        result = {
            'Kromozom': chrom,
            'Pozisyon': pos,
            'Referans': ref,
            'Alternatif': alt,
            'Okuma Derinliği': dp,
            'Varyant Yüzdesi': vaf,
            'Varyant Tipi': variant_type,
            'Varyant Etkisi': effect,
            'Mutasyon Tipi': mutation_type,
            'Protein Etkisi': protein_impact,
            'Protein Fonksiyonu': protein_function,
            'Klinik Etki': clinical_impact,
            'Gen': gene_info.get('Gen', '-'),
            'Kanser Geni': '+' if gene_info and gene_info.get('Kanser Geni') else '-',
            'Kanser Tipi': self.cancer_genes.get(gene_info.get('Gen', ''), '-'),
            'Transkript': gene_info.get('Transkript', '-'),
            'ClinVar': '+' if variant_annotations['ClinVar'].get('Bulundu', False) else '-',
            'ClinVar Önemi': variant_annotations['ClinVar'].get('Klinik Önemi', '-'),
            'ClinVar Durumu': variant_annotations['ClinVar'].get('İnceleme Durumu', '-'),
            'gnomAD Frekans': variant_annotations['gnomAD'].get('Allel Frekansı', '-'),
            'CADD Skoru': variant_annotations['MyVariant'].get('CADD Skoru', '-'),
            'SIFT': variant_annotations['MyVariant'].get('SIFT Tahmini', '-'),
            'PolyPhen': variant_annotations['MyVariant'].get('PolyPhen Tahmini', '-')
        }
        return result
class VCFAnalyzer(MultiDBVariantAnnotator):
    def __init__(self):
        super().__init__()
        self.root = tk.Tk()
        self.root.title("SDÜ NGS Panel Varyant Analizi")
        self.setup_gui()

    def setup_gui(self):
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Panel seçimi
        panel_frame = ttk.LabelFrame(self.main_frame, text="Panel Tipi", padding="5")
        panel_frame.pack(fill=tk.X, pady=5)
        
        self.panel_var = tk.StringVar(value="all")  # Varsayılan değer
        ttk.Radiobutton(panel_frame, text="Tüm Genler", 
                       variable=self.panel_var, value="all").pack(side=tk.LEFT)
        ttk.Radiobutton(panel_frame, text="Solid Tümör Paneli", 
                       variable=self.panel_var, value="solid").pack(side=tk.LEFT)
        ttk.Radiobutton(panel_frame, text="Akciğer Paneli", 
                       variable=self.panel_var, value="lung").pack(side=tk.LEFT)

        # Veritabanı seçimi
        db_frame = ttk.LabelFrame(self.main_frame, text="Veritabanı Sorguları (Opsiyonel)", padding="5")
        db_frame.pack(fill=tk.X, pady=5)
        
        self.clinvar_var = tk.BooleanVar(value=False)
        self.gnomad_var = tk.BooleanVar(value=False)
        self.myvariant_var = tk.BooleanVar(value=False)
        
        ttk.Checkbutton(db_frame, text="ClinVar", variable=self.clinvar_var).pack(side=tk.LEFT)
        ttk.Checkbutton(db_frame, text="gnomAD", variable=self.gnomad_var).pack(side=tk.LEFT)
        ttk.Checkbutton(db_frame, text="MyVariant.info", variable=self.myvariant_var).pack(side=tk.LEFT)
        
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

    def read_vcf(self, vcf_file):
        variants = []
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 10:  # FORMAT ve SAMPLE alanlarını kontrol et
                    format_fields = fields[8].split(':')
                    sample_fields = fields[9].split(':')
                    
                    # DP ve CLCAD2 indekslerini bul
                    dp_idx = format_fields.index('DP') if 'DP' in format_fields else -1
                    clcad2_idx = format_fields.index('CLCAD2') if 'CLCAD2' in format_fields else -1
                    
                    # Okuma derinliği (DP)
                    dp = '-'
                    if dp_idx != -1 and len(sample_fields) > dp_idx:
                        dp = sample_fields[dp_idx]

                    # Varyant yüzdesi (VAF)
                    vaf = '-'
                    if clcad2_idx != -1 and len(sample_fields) > clcad2_idx:
                        try:
                            ad_values = sample_fields[clcad2_idx].split(',')
                            if len(ad_values) >= 2:
                                ref_depth = int(ad_values[0])
                                alt_depth = int(ad_values[1])
                                total_depth = ref_depth + alt_depth
                                if total_depth > 0:
                                    vaf = f"{(alt_depth/total_depth)*100:.2f}%"
                        except:
                            pass

                    variants.append({
                        'CHROM': fields[0],
                        'POS': fields[1],
                        'ID': fields[2],
                        'REF': fields[3],
                        'ALT': fields[4],
                        'DP': dp,
                        'VAF': vaf
                    })
        return variants

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
                        variant['ID'],
                        variant['DP'],
                        variant['VAF']
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
            self.save_excel_report(df, patient_id, panel_type)
            self.create_visualizations(df, patient_id, panel_type)
            
            self.update_status(f"Analiz tamamlandı!", 100)
            
            completion_message = f"Toplam {len(df)} varyant analiz edildi."
            if panel_type != "all":
                panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                panel_variants = df[df['Gen'].isin(panel_genes)]
                completion_message += f"\nPanel genlerinde {len(panel_variants)} varyant bulundu."
            messagebox.showinfo("Analiz Tamamlandı", completion_message)
            
        except Exception as e:
            self.update_status(f"Hata: {str(e)}")
            messagebox.showerror("Hata", str(e))

    def save_excel_report(self, df, patient_id, panel_type):
        """Geliştirilmiş Excel raporu"""
        try:
            output_file = f"{patient_id}_varyant_analizi.xlsx"
            
            # Boş değerleri '-' ile değiştir
            df = df.fillna('-')
            
            with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
                # Excel formatları
                workbook = writer.book
                header_format = workbook.add_format({
                    'bold': True,
                    'bg_color': '#D7E4BC',
                    'border': 1
                })
                
                # Sütun sıralaması
                column_order = [
                    'Kromozom',
                    'Pozisyon',
                    'Referans',
                    'Alternatif',
                    'Okuma Derinliği',
                    'Varyant Yüzdesi',
                    'Varyant Tipi',
                    'Varyant Etkisi',
                    'Mutasyon Tipi',
                    'Protein Etkisi',
                    'Protein Fonksiyonu',
                    'Klinik Etki',
                    'Gen',
                    'Kanser Geni',
                    'Kanser Tipi',
                    'Transkript'
                ]
                
                # Veritabanı sorgularına göre ek sütunlar
                if self.clinvar_var.get():
                    column_order.extend(['ClinVar', 'ClinVar Önemi', 'ClinVar Durumu'])
                if self.gnomad_var.get():
                    column_order.append('gnomAD Frekans')
                if self.myvariant_var.get():
                    column_order.extend(['CADD Skoru', 'SIFT', 'PolyPhen'])
                
                # Tüm varyantlar
                df_ordered = df[column_order]
                df_ordered.to_excel(writer, sheet_name='Tüm Varyantlar', index=False)
                
                # Kanser genleri
                cancer_variants = df[df['Kanser Geni'] == '+'].sort_values(
                    by=['Klinik Etki', 'Gen'],
                    ascending=[False, True]
                )[column_order]
                if not cancer_variants.empty:
                    cancer_variants.to_excel(writer, sheet_name='Kanser Genleri', index=False)
                
                # Panel-spesifik varyantlar
                if panel_type != "all":
                    panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                    panel_variants = df[df['Gen'].isin(panel_genes)][column_order]
                    if not panel_variants.empty:
                        panel_variants.to_excel(writer, sheet_name='Panel Varyantları', index=False)
                
                # Yüksek etkili varyantlar
                high_impact = df[df['Klinik Etki'] == 'Yüksek'][column_order]
                if not high_impact.empty:
                    high_impact.to_excel(writer, sheet_name='Yüksek Etkili', index=False)
                
                # Özet istatistikler
                summary_data = {
                    'Metrik': [
                        'Toplam Varyant',
                        'Kanser Geni (+)',
                        'ClinVar (+)',
                        'Frameshift',
                        'SNV',
                        'İnsersiyon',
                        'Delesyon',
                        'Transition',
                        'Transversion',
                        'Yüksek Etki',
                        'Orta Etki',
                        'Düşük Etki'
                    ],
                    'Değer': [
                        len(df),
                        len(df[df['Kanser Geni'] == '+']),
                        len(df[df['ClinVar'] == '+']) if 'ClinVar' in df.columns else 0,
                        len(df[df['Varyant Etkisi'] == 'frameshift_variant']),
                        len(df[df['Varyant Tipi'] == 'SNV']),
                        len(df[df['Varyant Tipi'] == 'İnsersiyon']),
                        len(df[df['Varyant Tipi'] == 'Delesyon']),
                        len(df[df['Varyant Etkisi'] == 'transition']),
                        len(df[df['Varyant Etkisi'] == 'transversion']),
                        len(df[df['Klinik Etki'] == 'Yüksek']),
                        len(df[df['Klinik Etki'] == 'Orta']),
                        len(df[df['Klinik Etki'] == 'Düşük'])
                    ]
                }
                
                if panel_type != "all":
                    panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
                    panel_variants = df[df['Gen'].isin(panel_genes)]
                    summary_data['Metrik'].extend([
                        'Panel Genleri',
                        'Panel Oranı (%)'
                    ])
                    summary_data['Değer'].extend([
                        len(panel_variants),
                        round(len(panel_variants) / len(df) * 100, 2)
                    ])
                
                pd.DataFrame(summary_data).to_excel(writer, sheet_name='Özet', index=False)
                
                # Her sayfa için sütun genişliklerini ayarla
                for sheet_name in writer.sheets:
                    worksheet = writer.sheets[sheet_name]
                    for idx, col in enumerate(df_ordered.columns):
                        max_length = max(
                            df_ordered[col].astype(str).apply(len).max(),
                            len(col)
                        ) + 2
                        worksheet.set_column(idx, idx, min(max_length, 30))
            
            print(f"Excel raporu başarıyla kaydedildi: {output_file}")
            # Excel dosyasını otomatik aç
            if os.path.exists(output_file):
                os.startfile(output_file)
                
        except Exception as e:
            print(f"Excel raporu oluşturulurken hata oluştu: {str(e)}")
            print("DataFrame sütunları:", df.columns.tolist())
            messagebox.showerror("Hata", f"Excel raporu oluşturulamadı: {str(e)}")

    def create_visualizations(self, df, patient_id, panel_type):
        """Geliştirilmiş görselleştirmeler"""
        sns.set_theme(style="whitegrid")
        
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Varyant tipi dağılımı
        plt.subplot(2, 3, 1)
        sns.countplot(data=df, x='Varyant Tipi')
        plt.title('Varyant Tipi Dağılımı')
        plt.xticks(rotation=45)
        
        # 2. Gen dağılımı
        plt.subplot(2, 3, 2)
        if panel_type != "all":
            panel_genes = self.solid_panel_genes if panel_type == "solid" else self.lung_panel_genes
            panel_df = df[df['Gen'].isin(panel_genes)]
            if not panel_df.empty:
                sns.countplot(data=panel_df, y='Gen', 
                            order=panel_df['Gen'].value_counts().index)
                plt.title(f'{panel_type.upper()} Panel Genleri Dağılımı')
        else:
            cancer_df = df[df['Kanser Geni'] == '+']
            if not cancer_df.empty:
                sns.countplot(data=cancer_df, y='Gen',
                            order=cancer_df['Gen'].value_counts().head(15).index)
                plt.title('En Sık Görülen Kanser Genleri (Top 15)')
        
        # 3. Varyant etki dağılımı
        plt.subplot(2, 3, 3)
        sns.countplot(data=df, y='Varyant Etkisi')
        plt.title('Varyant Etki Dağılımı')
        
        # 4. Klinik etki dağılımı
        plt.subplot(2, 3, 4)
        sns.countplot(data=df, y='Klinik Etki')
        plt.title('Klinik Etki Dağılımı')
        
        # 5. Okuma derinliği ve VAF dağılımı
        plt.subplot(2, 3, 5)
        df['Okuma Derinliği'] = pd.to_numeric(df['Okuma Derinliği'].replace('-', '0'))
        df['Varyant Yüzdesi'] = pd.to_numeric(df['Varyant Yüzdesi'].str.rstrip('%').replace('-', '0'))
        plt.scatter(df['Okuma Derinliği'], df['Varyant Yüzdesi'], alpha=0.5)
        plt.xlabel('Okuma Derinliği')
        plt.ylabel('Varyant Yüzdesi (%)')
        plt.title('Okuma Derinliği vs Varyant Yüzdesi')
        
        plt.tight_layout()
        plt.savefig(f'{patient_id}_analiz_grafikleri.png', dpi=300, bbox_inches='tight')
        plt.close()

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
