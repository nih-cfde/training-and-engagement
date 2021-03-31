This is the JSON script for our demo workflow. It should be saved with a ".json" extension.

```
{
  "GWAS.sample_name": "coatColor",
  "GWAS.inputvcf": "pruned_coatColor_maf_geno.vcf",
  "GWAS.inputpheno": "coatColor.pheno",
  "GWAS.plink_missing_rates.inputminoralleles": "minor_alleles",
  "GWAS.plink_missing_rates.inped": "coatColor.ped",
  "GWAS.plink_missing_rates.inmap": "coatColor.map",
  "GWAS.plink_binary.inped": "coatColor.ped",
  "GWAS.plink_binary.inmap": "coatColor.map",
  "GWAS.plink_association.inputminoralleles": "minor_alleles",
  "GWAS.plink_association.inputbed": "coatColor.binary.bed",
  "GWAS.plink_association.inputfam": "coatColor.binary.fam",
  "GWAS.plink_association.inputbim": "coatColor.binary.bim",
  "GWAS.run_R.inadj": "coatColor.assoc.adjusted",
  "GWAS.run_R.inassoc": "coatColor.assoc"
}
```
