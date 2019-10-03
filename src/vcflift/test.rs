use super::*;
use std::fs::{self, File};

#[test]
fn test_header_lift() -> Result<(), LiftOverError> {
    let mut sample_header = &include_bytes!(
        "../../testfiles/annotations/gnomad.genomes.r2.1.1.sites.chr22.subset.header.txt"
    )[..];
    let reader = VCFReader::new(&mut sample_header)?;
    let vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new(),
    )?;

    let (_lifted_header, rewrite_target) = vcf_lift.lift_header(&reader.header)?;

    let expected = VCFHeaderRewriteTarget {
        info_alt: vec![
            b"number_a".to_vec(),
            b"AC".to_vec(),                //temporary
            b"AC_nfe_seu".to_vec(),        //temporary
            b"non_topmed_AC_amr".to_vec(), //temporary
            b"non_neuro_AC".to_vec(),      //temporary
        ]
        .into_iter()
        .collect(),
        info_ref: vec![b"number_r".to_vec()].into_iter().collect(),
        info_genotype: vec![b"number_g".to_vec()].into_iter().collect(),
        format_alt: vec![b"f_number_a".to_vec()].into_iter().collect(),
        format_ref: vec![b"f_number_r".to_vec()].into_iter().collect(),
        format_genotype: vec![b"f_number_g".to_vec()].into_iter().collect(),
        allele_count: vec![
            b"AC".to_vec(),
            b"AC_nfe_seu".to_vec(),
            b"non_topmed_AC_amr".to_vec(),
            b"non_neuro_AC".to_vec(),
        ]
        .into_iter()
        .collect(),
        allele_frequency: vec![
            b"AF".to_vec(),
            b"AF_nfe_seu".to_vec(),
            b"non_topmed_AF_amr".to_vec(),
            b"non_neuro_AF".to_vec(),
        ]
        .into_iter()
        .collect(),
        format_gt: false,
    };
    assert_eq!(rewrite_target, expected);

    Ok(())
}

#[test]
fn test_lift_vcf() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;
    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/lift-variant-test/lifttest/seq-a--to--seq-b.chain",
        "testfiles/lift-variant-test/sequence/seq-a.fa",
        "testfiles/lift-variant-test/sequence/seq-b.fa",
        VCFLiftOverParameters::new(),
    )?;
    let reader = File::open("testfiles/lift-variant-test/lifttest/seq-a.vcf")?;
    let success_writer = File::create("target/test-output/lifttest-success.vcf")?;
    let failed_writer = File::create("target/test-output/lifttest-fail.vcf")?;
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    //TODO: write test

    Ok(())
}

#[test]
fn test_lift_vcf2() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;
    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/lift-variant-test/lifttest/seq-a--to--seq-b.chain",
        "testfiles/lift-variant-test/sequence/seq-a.fa",
        "testfiles/lift-variant-test/sequence/seq-b.fa",
        VCFLiftOverParameters::new(),
    )?;
    let reader = File::open("testfiles/lift-variant-test/lifttest/seq-a.vcf")?;
    let success_writer = File::create("target/test-output/lifttest-success2.vcf")?;
    let failed_writer = File::create("target/test-output/lifttest-fail2.vcf")?;
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    //TODO: write test

    Ok(())
}

#[test]
fn test_lift_vcf_entry() -> Result<(), LiftOverError> {
    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new(),
    )?;

    let mut sample_header = &include_bytes!(
        "../../testfiles/annotations/gnomad.genomes.r2.0.2.sites.chr22.subset.header.all.vcf"
    )[..];
    let reader = VCFReader::new(&mut sample_header)?;

    let (_lifted_header, rewrite_target) = vcf_lift.lift_header(&reader.header)?;

    // rs533443276
    let original_entry = PartialVCFRecord::parse_vcf(1, b"22\t16581407\trs533443276\tT\tC\t3921.91\tPASS\tAC=13;AF=4.21914e-04;AN=30812;BaseQRankSum=-1.70200e+00;ClippingRankSum=-4.98000e-01;DP=679689;FS=5.21000e-01;InbreedingCoeff=-5.00000e-04;MQ=4.83700e+01;MQRankSum=-9.83000e-01;QD=8.89000e+00;ReadPosRankSum=2.05000e-01;SOR=7.54000e-01;VQSLOD=-7.01300e+00;VQSR_culprit=MQ;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|13;DP_HIST_ALT=0|0|1|1|2|3|0|3|1|0|2|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|0|0|1|0|0|0|3|1|4|2|1|1|0|0|0|0|0|0;GQ_HIST_ALL=5|14|12|30|29|29|57|108|103|265|355|288|1102|472|1059|626|1694|318|1626|7304;DP_HIST_ALL=16|60|142|576|1837|2704|3315|4005|1996|560|185|56|28|10|4|1|0|1|0|0;AB_HIST_ALL=0|0|0|0|1|0|0|0|3|1|4|2|1|1|0|0|0|0|0|0;AC_Male=5;AC_Female=8;AN_Male=17018;AN_Female=13794;AF_Male=2.93807e-04;AF_Female=5.79962e-04;GC_Male=8504,5,0;GC_Female=6889,8,0;GC_raw=15483,13,0;AC_raw=13;AN_raw=30992;GC=15393,13,0;AF_raw=4.19463e-04;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom=0;Hom_raw=0;AC_AFR=13;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_NFE=0;AC_OTH=0;AN_AFR=8692;AN_AMR=836;AN_ASJ=302;AN_EAS=1520;AN_FIN=3494;AN_NFE=14988;AN_OTH=980;AF_AFR=1.49563e-03;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;POPMAX=AFR;AC_POPMAX=13;AN_POPMAX=8692;AF_POPMAX=1.49563e-03;DP_MEDIAN=28;DREF_MEDIAN=1.25893e-32;GQ_MEDIAN=99;AB_MEDIAN=5.00000e-01;AS_RF=9.70323e-01;AS_FilterStatus=PASS;CSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs533443276|1||||SNV|1||||||||||||||||C:0.0004||C:0.0015|C:0|C:0|C:0|C:0||||||||||||||||||||||;GC_AFR=4333,13,0;GC_AMR=418,0,0;GC_ASJ=151,0,0;GC_EAS=760,0,0;GC_FIN=1747,0,0;GC_NFE=7494,0,0;GC_OTH=490,0,0;Hom_Male=0;Hom_Female=0\n")?;
    let mut expected_entry = PartialVCFRecord::parse_vcf(1, b"chr22\t15396556\trs533443276\tA\tG\t3921.91\tPASS\tAC=13;AF=4.21914e-04;AN=30812;BaseQRankSum=-1.70200e+00;ClippingRankSum=-4.98000e-01;DP=679689;FS=5.21000e-01;InbreedingCoeff=-5.00000e-04;MQ=4.83700e+01;MQRankSum=-9.83000e-01;QD=8.89000e+00;ReadPosRankSum=2.05000e-01;SOR=7.54000e-01;VQSLOD=-7.01300e+00;VQSR_culprit=MQ;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|13;DP_HIST_ALT=0|0|1|1|2|3|0|3|1|0|2|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|0|0|1|0|0|0|3|1|4|2|1|1|0|0|0|0|0|0;GQ_HIST_ALL=5|14|12|30|29|29|57|108|103|265|355|288|1102|472|1059|626|1694|318|1626|7304;DP_HIST_ALL=16|60|142|576|1837|2704|3315|4005|1996|560|185|56|28|10|4|1|0|1|0|0;AB_HIST_ALL=0|0|0|0|1|0|0|0|3|1|4|2|1|1|0|0|0|0|0|0;AC_Male=5;AC_Female=8;AN_Male=17018;AN_Female=13794;AF_Male=2.93807e-04;AF_Female=5.79962e-04;GC_Male=8504,5,0;GC_Female=6889,8,0;GC_raw=15483,13,0;AC_raw=13;AN_raw=30992;GC=15393,13,0;AF_raw=4.19463e-04;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom=0;Hom_raw=0;AC_AFR=13;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_NFE=0;AC_OTH=0;AN_AFR=8692;AN_AMR=836;AN_ASJ=302;AN_EAS=1520;AN_FIN=3494;AN_NFE=14988;AN_OTH=980;AF_AFR=1.49563e-03;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;POPMAX=AFR;AC_POPMAX=13;AN_POPMAX=8692;AF_POPMAX=1.49563e-03;DP_MEDIAN=28;DREF_MEDIAN=1.25893e-32;GQ_MEDIAN=99;AB_MEDIAN=5.00000e-01;AS_RF=9.70323e-01;AS_FilterStatus=PASS;CSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs533443276|1||||SNV|1||||||||||||||||C:0.0004||C:0.0015|C:0|C:0|C:0|C:0||||||||||||||||||||||;GC_AFR=4333,13,0;GC_AMR=418,0,0;GC_ASJ=151,0,0;GC_EAS=760,0,0;GC_FIN=1747,0,0;GC_NFE=7494,0,0;GC_OTH=490,0,0;Hom_Male=0;Hom_Female=0;ORIGINAL_CHROM=22;ORIGINAL_POS=16581407;ORIGINAL_STRAND=-\n")?;
    expected_entry.raw = original_entry.raw;
    expected_entry.original_unparsed_info = original_entry.original_unparsed_info;
    let lifted_entry = vcf_lift.lift_record(&original_entry, &rewrite_target)?;
    assert_eq!(
        lifted_entry,
        VCFLiftOverResult::Succeeded(vec![VCFRecordWrapper::Partial(expected_entry)])
    );

    // rs387906933
    let original_entry = PartialVCFRecord::parse_vcf(1, b"22\t51137225\trs387906933\tC\tG\t209.47\tPASS\tAC=1;AF=3.23060e-05;AN=30954;BaseQRankSum=-1.38000e+00;ClippingRankSum=-8.73000e-01;DP=673396;FS=7.14500e+00;InbreedingCoeff=2.40000e-03;MQ=6.00000e+01;MQRankSum=-6.52000e-01;QD=6.98000e+00;ReadPosRankSum=1.02000e+00;SOR=1.41900e+00;VQSLOD=1.30000e+01;VQSR_culprit=MQ;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST_ALL=2|2|2|7|8|20|50|119|106|296|431|314|1353|417|1265|629|1907|267|1758|6542;DP_HIST_ALL=4|8|172|759|2308|3127|3756|3789|1086|298|105|48|18|7|1|3|0|0|2|1;AB_HIST_ALL=0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0;AC_Male=1;AC_Female=0;AN_Male=17106;AN_Female=13848;AF_Male=5.84590e-05;AF_Female=0.00000e+00;GC_Male=8552,1,0;GC_Female=6924,0,0;GC_raw=15494,1,0;AC_raw=1;AN_raw=30990;GC=15476,1,0;AF_raw=3.22685e-05;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom=0;Hom_raw=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=1;AC_FIN=0;AC_NFE=0;AC_OTH=0;AN_AFR=8728;AN_AMR=836;AN_ASJ=302;AN_EAS=1622;AN_FIN=3494;AN_NFE=14992;AN_OTH=980;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=6.16523e-04;AF_FIN=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;POPMAX=EAS;AC_POPMAX=1;AN_POPMAX=1622;AF_POPMAX=6.16523e-04;DP_MEDIAN=30;DREF_MEDIAN=6.30957e-28;GQ_MEDIAN=99;AB_MEDIAN=3.66667e-01;AS_RF=8.98417e-01;AS_FilterStatus=PASS;CSQ=G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000262795|protein_coding|13/23||ENST00000262795.3:c.1654C>G|ENSP00000442518.1:p.Arg552Gly|1654|1654|552|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|YES||||ENSP00000442518||F8TCV3&F2Z3L0|UPI0000DD85FB|1||possibly_damaging(0.586)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||,G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000414786|protein_coding|12/23||ENST00000414786.2:c.1564C>G|ENSP00000464552.1:p.Arg522Gly|1791|1564|522|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|||||ENSP00000464552||M0QWZ9&F8TCV3|UPI000268B48A|1||possibly_damaging(0.722)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||,G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000445220|protein_coding|12/23||ENST00000445220.2:c.1609C>G|ENSP00000446078.1:p.Arg537Gly|1609|1609|537|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|||||ENSP00000446078||F8TCV3&F5GZA1|UPI000206575B|1|deleterious(0)|possibly_damaging(0.468)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||;GC_AFR=4364,0,0;GC_AMR=418,0,0;GC_ASJ=151,0,0;GC_EAS=810,1,0;GC_FIN=1747,0,0;GC_NFE=7496,0,0;GC_OTH=490,0,0;Hom_Male=0;Hom_Female=0\n")?;
    let mut expected_entry = PartialVCFRecord::parse_vcf(1, b"chr22\t50698797\trs387906933\tC\tG\t209.47\tPASS\tAC=1;AF=3.23060e-05;AN=30954;BaseQRankSum=-1.38000e+00;ClippingRankSum=-8.73000e-01;DP=673396;FS=7.14500e+00;InbreedingCoeff=2.40000e-03;MQ=6.00000e+01;MQRankSum=-6.52000e-01;QD=6.98000e+00;ReadPosRankSum=1.02000e+00;SOR=1.41900e+00;VQSLOD=1.30000e+01;VQSR_culprit=MQ;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST_ALL=2|2|2|7|8|20|50|119|106|296|431|314|1353|417|1265|629|1907|267|1758|6542;DP_HIST_ALL=4|8|172|759|2308|3127|3756|3789|1086|298|105|48|18|7|1|3|0|0|2|1;AB_HIST_ALL=0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0;AC_Male=1;AC_Female=0;AN_Male=17106;AN_Female=13848;AF_Male=5.84590e-05;AF_Female=0.00000e+00;GC_Male=8552,1,0;GC_Female=6924,0,0;GC_raw=15494,1,0;AC_raw=1;AN_raw=30990;GC=15476,1,0;AF_raw=3.22685e-05;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom=0;Hom_raw=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=1;AC_FIN=0;AC_NFE=0;AC_OTH=0;AN_AFR=8728;AN_AMR=836;AN_ASJ=302;AN_EAS=1622;AN_FIN=3494;AN_NFE=14992;AN_OTH=980;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=6.16523e-04;AF_FIN=0.00000e+00;AF_NFE=0.00000e+00;AF_OTH=0.00000e+00;POPMAX=EAS;AC_POPMAX=1;AN_POPMAX=1622;AF_POPMAX=6.16523e-04;DP_MEDIAN=30;DREF_MEDIAN=6.30957e-28;GQ_MEDIAN=99;AB_MEDIAN=3.66667e-01;AS_RF=8.98417e-01;AS_FilterStatus=PASS;CSQ=G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000262795|protein_coding|13/23||ENST00000262795.3:c.1654C>G|ENSP00000442518.1:p.Arg552Gly|1654|1654|552|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|YES||||ENSP00000442518||F8TCV3&F2Z3L0|UPI0000DD85FB|1||possibly_damaging(0.586)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||,G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000414786|protein_coding|12/23||ENST00000414786.2:c.1564C>G|ENSP00000464552.1:p.Arg522Gly|1791|1564|522|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|||||ENSP00000464552||M0QWZ9&F8TCV3|UPI000268B48A|1||possibly_damaging(0.722)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||,G|missense_variant|MODERATE|SHANK3|ENSG00000251322|Transcript|ENST00000445220|protein_coding|12/23||ENST00000445220.2:c.1609C>G|ENSP00000446078.1:p.Arg537Gly|1609|1609|537|R/G|Cgg/Ggg|rs387906933|1||1||SNV|1|HGNC|14294|||||ENSP00000446078||F8TCV3&F5GZA1|UPI000206575B|1|deleterious(0)|possibly_damaging(0.468)|hmmpanther:PTHR24135&hmmpanther:PTHR24135:SF4|||||||||||||||||||pathogenic||1|||||||||;GC_AFR=4364,0,0;GC_AMR=418,0,0;GC_ASJ=151,0,0;GC_EAS=810,1,0;GC_FIN=1747,0,0;GC_NFE=7496,0,0;GC_OTH=490,0,0;Hom_Male=0;Hom_Female=0;ORIGINAL_CHROM=22;ORIGINAL_POS=51137225;ORIGINAL_STRAND=+\n")?;
    expected_entry.raw = original_entry.raw;
    expected_entry.original_unparsed_info = original_entry.original_unparsed_info;
    let lifted_entry = vcf_lift.lift_record(&original_entry, &rewrite_target)?;
    assert_eq!(
        lifted_entry,
        VCFLiftOverResult::Succeeded(vec![VCFRecordWrapper::Partial(expected_entry)])
    );

    // rs267608319
    let original_entry = PartialVCFRecord::parse_vcf(1, b"22\t42522751\t.\tC\tT\t.\tPASS\t.\n")?;
    let mut expected_entry = PartialVCFRecord::parse_vcf(1, b"chr22\t42126749\t.\tC\tT\t.\tPASS\tORIGINAL_CHROM=22;ORIGINAL_POS=42522751;ORIGINAL_STRAND=+\n")?;
    expected_entry.raw = original_entry.raw;
    expected_entry.original_unparsed_info = original_entry.original_unparsed_info;
    let lifted_entry = vcf_lift.lift_record(&original_entry, &rewrite_target)?;
    assert_eq!(
        lifted_entry,
        VCFLiftOverResult::Succeeded(vec![VCFRecordWrapper::Partial(expected_entry)])
    );

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new().do_not_prefer_cis_contig_when_multimap(true),
    )?;
    let original_entry = PartialVCFRecord::parse_vcf(1, b"22\t42522751\t.\tC\tT\t.\tPASS\t.\n")?;
    let mut expected_entry = PartialVCFRecord::parse_vcf(
        1,
        b"22\t42522751\t.\tC\tT\t.\tPASS\tFAILED_REASON=MULTIMAP;MULTIMAP=2\n",
    )?;
    expected_entry.raw = original_entry.raw;
    expected_entry.original_unparsed_info = original_entry.original_unparsed_info;
    let lifted_entry = vcf_lift.lift_record(&original_entry, &rewrite_target)?;
    assert_eq!(
        lifted_entry,
        VCFLiftOverResult::Failed(Box::new(expected_entry))
    );

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new()
            .do_not_prefer_cis_contig_when_multimap(true)
            .allow_multimap(true),
    )?;
    let original_entry = PartialVCFRecord::parse_vcf(1, b"22\t42522751\t.\tC\tT\t.\tPASS\t.\n")?;
    let mut expected_entry1 =
        PartialVCFRecord::parse_vcf(1, b"chr22\t42126749\t.\tC\tT\t.\tPASS\tORIGINAL_CHROM=22;ORIGINAL_POS=42522751;ORIGINAL_STRAND=+;MULTIMAP=2\n")?;
    expected_entry1.raw = original_entry.raw;
    expected_entry1.original_unparsed_info = original_entry.original_unparsed_info;
    let mut expected_entry2 = PartialVCFRecord::parse_vcf(
        1,
        b"chr22_KI270928v1_alt\t49090\t.\tC\tT\t.\tPASS\tORIGINAL_CHROM=22;ORIGINAL_POS=42522751;ORIGINAL_STRAND=+;MULTIMAP=2\n",
    )?;
    expected_entry2.raw = original_entry.raw;
    expected_entry2.original_unparsed_info = original_entry.original_unparsed_info;

    let lifted_entry = vcf_lift.lift_record(&original_entry, &rewrite_target)?;
    if let VCFLiftOverResult::Succeeded(mut succeeded_records) = lifted_entry {
        succeeded_records.sort_by_key(|x| x.position());
        assert_eq!(
            succeeded_records,
            vec![
                VCFRecordWrapper::Partial(expected_entry2),
                VCFRecordWrapper::Partial(expected_entry1)
            ]
        );
    } else {
        panic!();
    }

    Ok(())
}

#[test]
fn test_lift_1000genomes() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new(),
    )?;
    let reader = flate2::read::MultiGzDecoder::new(File::open("testfiles/annotations/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.subset.vcf.gz")?);
    let success_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-1000genomes-success.vcf.gz")?,
        flate2::Compression::default(),
    );
    let failed_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-1000genomes-fail.vcf.gz")?,
        flate2::Compression::default(),
    );
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    //TODO: write test

    Ok(())
}

#[test]
fn test_real_data_gnomd_2_0_2() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new(),
    )?;
    let reader = flate2::read::MultiGzDecoder::new(File::open(
        "testfiles/annotations/gnomad.genomes.r2.0.2.sites.chr22.subset.vcf.bgz",
    )?);
    let success_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-gnomad-2.0.2-success.vcf.gz")?,
        flate2::Compression::default(),
    );
    let failed_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-gnomad-2.0.2-fail.vcf.gz")?,
        flate2::Compression::default(),
    );
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    Ok(())
}

#[test]
fn test_real_data_gnomd_2_1_1() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new().allow_multimap(true),
    )?;
    let reader = flate2::read::MultiGzDecoder::new(File::open(
        "testfiles/annotations/gnomad.genomes.r2.1.1.sites.chr22.subset.merged.vcf.gz",
    )?);
    let success_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-gnomad-2.1.1-success.vcf.gz")?,
        flate2::Compression::default(),
    );
    let failed_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-gnomad-2.1.1-fail.vcf.gz")?,
        flate2::Compression::default(),
    );
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    Ok(())
}

#[test]
fn test_real_data_clinvar() -> Result<(), LiftOverError> {
    fs::create_dir_all("target/test-output/liftvcf")?;

    let mut vcf_lift = VCFLiftOver::load(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
        "testfiles/genome/hg19-chr22.fa",
        "testfiles/genome/hg38-chr22.fa",
        VCFLiftOverParameters::new().do_not_swap_ref_alt(true),
    )?;
    let reader = flate2::read::MultiGzDecoder::new(File::open(
        "testfiles/annotations/clinvar_20190722.chr22.vcf.gz",
    )?);
    let success_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-clinvar_20190722.chr22-success.vcf.gz")?,
        flate2::Compression::default(),
    );
    let failed_writer = flate2::write::GzEncoder::new(
        File::create("target/test-output/liftvcf/lifttest-clinvar_20190722.chr22-fail.vcf.gz")?,
        flate2::Compression::default(),
    );
    vcf_lift.lift_vcf(reader, success_writer, failed_writer)?;

    Ok(())
}

#[test]
fn test_swap_ref_alt() -> Result<(), LiftOverError> {
    let rewrite_target = VCFHeaderRewriteTarget {
        info_ref: vec![b"I_R".to_vec()].into_iter().collect(),
        info_alt: vec![b"I_A".to_vec()].into_iter().collect(),
        info_genotype: vec![b"I_G".to_vec()].into_iter().collect(),
        allele_frequency: vec![b"AF".to_vec(), b"AF_EAS".to_vec()]
            .into_iter()
            .collect(),
        allele_count: vec![b"AC".to_vec(), b"AC_EAS".to_vec()]
            .into_iter()
            .collect(),
        format_ref: vec![b"F_R".to_vec()].into_iter().collect(),
        format_alt: vec![b"F_A".to_vec()].into_iter().collect(),
        format_genotype: vec![b"F_G".to_vec()].into_iter().collect(),
        format_gt: true,
    };
    let param = VCFLiftOverParameters::new();

    let original_record = PartialVCFRecord::parse_vcf(0, b"22\t17590404\t340611\tG\tA,T\t.\t.\tI_R=G1,A1,T1;I_A=A2,T2\tGT:F_R:F_A\t0/1:G3,A3,T3:A4,T4\t1|2:G5,A5,T5:A6,T6\n")?;
    let lifted_variant = LiftedVariant {
        chromosome: "20".to_string(),
        position: 1999,
        strand: crate::chain::Strand::Reverse,
        original_reference: b"G".to_vec(),
        reference: b"G".to_vec(),
        alternative: vec![b"A".to_vec(), b"T".to_vec()],
        reference_changed: false,
    };

    let expected = PartialVCFRecord::parse_vcf(0, b"20\t2000\t340611\tG\tA,T\t.\t.\tI_R=G1,A1,T1;I_A=A2,T2;ORIGINAL_CHROM=22;ORIGINAL_POS=17590404;ORIGINAL_STRAND=-\tGT:F_R:F_A\t0/1:G3,A3,T3:A4,T4\t1|2:G5,A5,T5:A6,T6\n")?.complete_parse()?;
    assert_eq!(
        swap_ref_alt(&lifted_variant, &original_record, &param, &rewrite_target)?,
        expected
    );

    let lifted_variant = LiftedVariant {
        chromosome: "10".to_string(),
        position: 1000,
        strand: crate::chain::Strand::Forward,
        original_reference: b"G".to_vec(),
        reference: b"A".to_vec(),
        alternative: vec![b"A".to_vec(), b"T".to_vec()],
        reference_changed: true,
    };
    let expected = PartialVCFRecord::parse_vcf(0, b"10\t1001\t340611\tA\tT,G\t.\t.\tI_R=A1,T1,G1;I_A=T2,.;ORIGINAL_CHROM=22;ORIGINAL_POS=17590404;ORIGINAL_STRAND=+;REF_CHANGED;ORIGINAL_REF=G\tGT:F_R:F_A\t2/0:A3,T3,G3:T4,.\t0|1:A5,T5,G5:T6,.\n")?.complete_parse()?;
    let swapped = swap_ref_alt(&lifted_variant, &original_record, &param, &rewrite_target)?;
    assert_eq!(swapped, expected);
    assert_eq!(
        merge_to_vcf(&lifted_variant, &original_record, &param, &rewrite_target)?,
        VCFRecordWrapper::Complete(expected)
    );

    let lifted_variant = LiftedVariant {
        chromosome: "15".to_string(),
        position: 9,
        strand: crate::chain::Strand::Forward,
        original_reference: b"G".to_vec(),
        reference: b"C".to_vec(),
        alternative: vec![b"A".to_vec(), b"T".to_vec()],
        reference_changed: true,
    };
    let expected = PartialVCFRecord::parse_vcf(0, b"15\t10\t340611\tC\tA,T,G\t.\t.\tI_R=.,A1,T1,G1;I_A=A2,T2,.;ORIGINAL_CHROM=22;ORIGINAL_POS=17590404;ORIGINAL_STRAND=+;REF_CHANGED;ORIGINAL_REF=G\tGT:F_R:F_A\t3/1:.,A3,T3,G3:A4,T4,.\t1|2:.,A5,T5,G5:A6,T6,.\n")?.complete_parse()?;
    let swapped = swap_ref_alt(&lifted_variant, &original_record, &param, &rewrite_target)?;
    assert_eq!(swapped, expected);
    assert_eq!(
        merge_to_vcf(&lifted_variant, &original_record, &param, &rewrite_target)?,
        VCFRecordWrapper::Complete(expected)
    );
    Ok(())
}
