use super::*;

#[allow(clippy::unreadable_literal)]
#[test]
fn test_parse_partial_vcf() {
    // samples
    let data = b"22\t16914962\trs9306245\tC\tA,T\t100\tPASS\tAC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC\tGT:X:Y\t0|2:.:1,2\t0|0:1,2,3:1,2\t0|2:.:.\t0|2:1,2:1,2,3\t0|2:x,y:x,y,z,w\n";
    let record = PartialVCFRecord::parse_vcf(1, data).unwrap();
    assert_eq!(
        record,
        PartialVCFRecord {
            raw: data,
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 16914962,
            id: Cow::Borrowed(b"rs9306245"),
            reference: Cow::Borrowed(b"C"),
            alternative: vec![Cow::Borrowed(b"A"), Cow::Borrowed(b"T")],
            qual: Cow::Borrowed(b"100"),
            filter: Cow::Borrowed(b"PASS"),
            unparsed_info: Cow::Borrowed(b"AC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC"),
            original_unparsed_info: b"AC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC",
            other: b"GT:X:Y\t0|2:.:1,2\t0|0:1,2,3:1,2\t0|2:.:.\t0|2:1,2:1,2,3\t0|2:x,y:x,y,z,w",
        }
    );
    let mut write_data = Vec::new();
    record.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);
    let complete_data = record.complete_parse().unwrap();
    assert_eq!(
        complete_data,
        CompleteVCFRecord {
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 16914962,
            id: Cow::Borrowed(b"rs9306245"),
            reference: Cow::Borrowed(b"C"),
            alternative: vec![Cow::Borrowed(b"A"), Cow::Borrowed(b"T")],
            qual: Cow::Borrowed(b"100"),
            filter: Cow::Borrowed(b"PASS"),
            info: vec![
                (
                    Cow::Borrowed(b"AC"),
                    vec![Cow::Borrowed(b"161"), Cow::Borrowed(b"1659")]
                ),
                (
                    Cow::Borrowed(b"AF"),
                    vec![Cow::Borrowed(b"0.0321486"), Cow::Borrowed(b"0.33127")]
                ),
                (Cow::Borrowed(b"AN"), vec![Cow::Borrowed(b"5008")]),
                (Cow::Borrowed(b"MULTI_ALLELIC"), vec![])
            ],
            format: vec![
                Cow::Borrowed(b"GT"),
                Cow::Borrowed(b"X"),
                Cow::Borrowed(b"Y")
            ],
            call: vec![
                vec![
                    vec![Cow::Borrowed(b"0|2")],
                    vec![Cow::Borrowed(b".")],
                    vec![Cow::Borrowed(b"1"), Cow::Borrowed(b"2")]
                ],
                vec![
                    vec![Cow::Borrowed(b"0|0")],
                    vec![
                        Cow::Borrowed(b"1"),
                        Cow::Borrowed(b"2"),
                        Cow::Borrowed(b"3")
                    ],
                    vec![Cow::Borrowed(b"1"), Cow::Borrowed(b"2")]
                ],
                vec![
                    vec![Cow::Borrowed(b"0|2")],
                    vec![Cow::Borrowed(b".")],
                    vec![Cow::Borrowed(b".")]
                ],
                vec![
                    vec![Cow::Borrowed(b"0|2")],
                    vec![Cow::Borrowed(b"1"), Cow::Borrowed(b"2")],
                    vec![
                        Cow::Borrowed(b"1"),
                        Cow::Borrowed(b"2"),
                        Cow::Borrowed(b"3")
                    ]
                ],
                vec![
                    vec![Cow::Borrowed(b"0|2")],
                    vec![Cow::Borrowed(b"x"), Cow::Borrowed(b"y")],
                    vec![
                        Cow::Borrowed(b"x"),
                        Cow::Borrowed(b"y"),
                        Cow::Borrowed(b"z"),
                        Cow::Borrowed(b"w")
                    ]
                ]
            ],
        }
    );
    let mut write_data = Vec::new();
    complete_data.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);

    // clinVar
    let data = b"22\t17590404\t340611\tG\tA\t.\t.\tAF_ESP=0.02504;AF_EXAC=0.01002;AF_TGP=0.02196;ALLELEID=352044;CLNDISDB=MedGen:C4310803,OMIM:613953|MedGen:CN239217\n";
    let record = PartialVCFRecord::parse_vcf(1, data).unwrap();

    assert_eq!(record, PartialVCFRecord{
            raw: data,
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 17590404,
            id: Cow::Borrowed(b"340611"),
            reference: Cow::Borrowed(b"G"),
            alternative: vec![Cow::Borrowed(b"A")],
            qual: Cow::Borrowed(b"."),
            filter: Cow::Borrowed(b"."),
            unparsed_info: Cow::Borrowed(b"AF_ESP=0.02504;AF_EXAC=0.01002;AF_TGP=0.02196;ALLELEID=352044;CLNDISDB=MedGen:C4310803,OMIM:613953|MedGen:CN239217"),
            original_unparsed_info: b"AF_ESP=0.02504;AF_EXAC=0.01002;AF_TGP=0.02196;ALLELEID=352044;CLNDISDB=MedGen:C4310803,OMIM:613953|MedGen:CN239217",
            other: b"",
        });
    let mut write_data = Vec::new();
    record.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);
    let complete_data = record.complete_parse().unwrap();
    assert_eq!(
        complete_data,
        CompleteVCFRecord {
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 17590404,
            id: Cow::Borrowed(b"340611"),
            reference: Cow::Borrowed(b"G"),
            alternative: vec![Cow::Borrowed(b"A")],
            qual: Cow::Borrowed(b"."),
            filter: Cow::Borrowed(b"."),
            info: vec![
                (Cow::Borrowed(b"AF_ESP"), vec![Cow::Borrowed(b"0.02504")]),
                (Cow::Borrowed(b"AF_EXAC"), vec![Cow::Borrowed(b"0.01002")]),
                (Cow::Borrowed(b"AF_TGP"), vec![Cow::Borrowed(b"0.02196")]),
                (Cow::Borrowed(b"ALLELEID"), vec![Cow::Borrowed(b"352044")]),
                (
                    Cow::Borrowed(b"CLNDISDB"),
                    vec![
                        Cow::Borrowed(b"MedGen:C4310803"),
                        Cow::Borrowed(b"OMIM:613953|MedGen:CN239217")
                    ]
                )
            ],
            format: vec![],
            call: vec![],
        }
    );
    let mut write_data = Vec::new();
    complete_data.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);

    // gnomAD 2.0.2
    let data = b"22\t16056854\trs58796206\tG\tGA,A\t9925077.10\tPASS\tAC=12072,1;AF=4.03611e-01,3.34336e-05;AN=29910;BaseQRankSum=-2.31000e-01;ClippingRankSum=-6.00000e-03;DP=710017;FS=6.12600e+00;InbreedingCoeff=-9.36000e-02;MQ=5.09600e+01;MQRankSum=-4.28700e+00;QD=2.00400e+01;ReadPosRankSum=-4.10000e-01;SOR=1.05100e+00;VQSLOD=-6.57000e+00;VQSR_culprit=MQRankSum;GQ_HIST_ALT=15|19|12|23|42|29|53|93|59|124|167|89|174|172|98|178|139|82|116|8511,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=6|73|421|1185|1785|1676|1300|844|663|514|428|347|263|214|134|120|69|49|35|17,0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|4|9|41|145|387|735|997|1373|1016|1051|570|496|440|365|258|109|31|17|0,0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0;GQ_HIST_ALL=173|42|29|68|96|71|158|218|162|337|412|281|699|353|549|396|680|182|574|9920;DP_HIST_ALL=12|102|557|1623|2757|2714|2369|1963|979|567|442|356|265|216|135|120|69|49|35|18;AB_HIST_ALL=0|4|10|41|145|387|735|995|1373|1018|1052|570|498|440|364|259|109|31|17|0;AC_Male=6665,1;AC_Female=5407,0;AN_Male=16504;AN_Female=13406;AF_Male=4.03841e-01,6.05914e-05;AF_Female=4.03327e-01,0.00000e+00;GC_Male=2756,4321,1172,1,0,0;GC_Female=2189,3619,894,0,0,0;GC_raw=5200,8044,2151,1,0,0;AC_raw=12346,1;AN_raw=30800;GC=4945,7940,2066,1,0,0;AF_raw=4.00844e-01,3.24675e-05;Hom_AFR=631,0;Hom_AMR=106,0;Hom_ASJ=20,0;Hom_EAS=304,0;Hom_FIN=193,0;Hom_NFE=767,0;Hom_OTH=45,0;Hom=2066,0;Hom_raw=2151,0;AC_AFR=3492,0;AC_AMR=442,0;AC_ASJ=91,0;AC_EAS=1019,0;AC_FIN=1460,0;AC_NFE=5207,1;AC_OTH=361,0;AN_AFR=8494;AN_AMR=804;AN_ASJ=260;AN_EAS=1590;AN_FIN=3484;AN_NFE=14330;AN_OTH=948;AF_AFR=4.11114e-01,0.00000e+00;AF_AMR=5.49751e-01,0.00000e+00;AF_ASJ=3.50000e-01,0.00000e+00;AF_EAS=6.40881e-01,0.00000e+00;AF_FIN=4.19059e-01,0.00000e+00;AF_NFE=3.63364e-01,6.97837e-05;AF_OTH=3.80802e-01,0.00000e+00;STAR_AC=4;STAR_AC_raw=5;STAR_Hom=1;POPMAX=EAS,NFE;AC_POPMAX=1019,1;AN_POPMAX=1590,14330;AF_POPMAX=6.40881e-01,6.97837e-05;DP_MEDIAN=30,30;DREF_MEDIAN=5.01187e-55,2.51189e-46;GQ_MEDIAN=99,99;AB_MEDIAN=4.61538e-01,6.00000e-01;AS_RF=7.75613e-01,7.57709e-01;AS_FilterStatus=PASS,PASS;CSQ=A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||1||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|intergenic_variant|MODIFIER||||||||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||;GC_AFR=1386,2230,631,0,0,0;GC_AMR=66,230,106,0,0,0;GC_ASJ=59,51,20,0,0,0;GC_EAS=80,411,304,0,0,0;GC_FIN=475,1074,193,0,0,0;GC_NFE=2721,3673,767,1,0,0;GC_OTH=158,271,45,0,0,0;Hom_Male=1172,0;Hom_Female=894,0\n";
    let record = PartialVCFRecord::parse_vcf(1, data).unwrap();
    assert_eq!(record, PartialVCFRecord{
            raw: data,
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 16056854,
            id: Cow::Borrowed(b"rs58796206"),
            reference: Cow::Borrowed(b"G"),
            alternative: vec![Cow::Borrowed(b"GA"), Cow::Borrowed(b"A")],
            qual: Cow::Borrowed(b"9925077.10"),
            filter: Cow::Borrowed(b"PASS"),
            unparsed_info: Cow::Borrowed(b"AC=12072,1;AF=4.03611e-01,3.34336e-05;AN=29910;BaseQRankSum=-2.31000e-01;ClippingRankSum=-6.00000e-03;DP=710017;FS=6.12600e+00;InbreedingCoeff=-9.36000e-02;MQ=5.09600e+01;MQRankSum=-4.28700e+00;QD=2.00400e+01;ReadPosRankSum=-4.10000e-01;SOR=1.05100e+00;VQSLOD=-6.57000e+00;VQSR_culprit=MQRankSum;GQ_HIST_ALT=15|19|12|23|42|29|53|93|59|124|167|89|174|172|98|178|139|82|116|8511,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=6|73|421|1185|1785|1676|1300|844|663|514|428|347|263|214|134|120|69|49|35|17,0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|4|9|41|145|387|735|997|1373|1016|1051|570|496|440|365|258|109|31|17|0,0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0;GQ_HIST_ALL=173|42|29|68|96|71|158|218|162|337|412|281|699|353|549|396|680|182|574|9920;DP_HIST_ALL=12|102|557|1623|2757|2714|2369|1963|979|567|442|356|265|216|135|120|69|49|35|18;AB_HIST_ALL=0|4|10|41|145|387|735|995|1373|1018|1052|570|498|440|364|259|109|31|17|0;AC_Male=6665,1;AC_Female=5407,0;AN_Male=16504;AN_Female=13406;AF_Male=4.03841e-01,6.05914e-05;AF_Female=4.03327e-01,0.00000e+00;GC_Male=2756,4321,1172,1,0,0;GC_Female=2189,3619,894,0,0,0;GC_raw=5200,8044,2151,1,0,0;AC_raw=12346,1;AN_raw=30800;GC=4945,7940,2066,1,0,0;AF_raw=4.00844e-01,3.24675e-05;Hom_AFR=631,0;Hom_AMR=106,0;Hom_ASJ=20,0;Hom_EAS=304,0;Hom_FIN=193,0;Hom_NFE=767,0;Hom_OTH=45,0;Hom=2066,0;Hom_raw=2151,0;AC_AFR=3492,0;AC_AMR=442,0;AC_ASJ=91,0;AC_EAS=1019,0;AC_FIN=1460,0;AC_NFE=5207,1;AC_OTH=361,0;AN_AFR=8494;AN_AMR=804;AN_ASJ=260;AN_EAS=1590;AN_FIN=3484;AN_NFE=14330;AN_OTH=948;AF_AFR=4.11114e-01,0.00000e+00;AF_AMR=5.49751e-01,0.00000e+00;AF_ASJ=3.50000e-01,0.00000e+00;AF_EAS=6.40881e-01,0.00000e+00;AF_FIN=4.19059e-01,0.00000e+00;AF_NFE=3.63364e-01,6.97837e-05;AF_OTH=3.80802e-01,0.00000e+00;STAR_AC=4;STAR_AC_raw=5;STAR_Hom=1;POPMAX=EAS,NFE;AC_POPMAX=1019,1;AN_POPMAX=1590,14330;AF_POPMAX=6.40881e-01,6.97837e-05;DP_MEDIAN=30,30;DREF_MEDIAN=5.01187e-55,2.51189e-46;GQ_MEDIAN=99,99;AB_MEDIAN=4.61538e-01,6.00000e-01;AS_RF=7.75613e-01,7.57709e-01;AS_FilterStatus=PASS,PASS;CSQ=A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||1||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|intergenic_variant|MODIFIER||||||||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||;GC_AFR=1386,2230,631,0,0,0;GC_AMR=66,230,106,0,0,0;GC_ASJ=59,51,20,0,0,0;GC_EAS=80,411,304,0,0,0;GC_FIN=475,1074,193,0,0,0;GC_NFE=2721,3673,767,1,0,0;GC_OTH=158,271,45,0,0,0;Hom_Male=1172,0;Hom_Female=894,0"),
            original_unparsed_info: b"AC=12072,1;AF=4.03611e-01,3.34336e-05;AN=29910;BaseQRankSum=-2.31000e-01;ClippingRankSum=-6.00000e-03;DP=710017;FS=6.12600e+00;InbreedingCoeff=-9.36000e-02;MQ=5.09600e+01;MQRankSum=-4.28700e+00;QD=2.00400e+01;ReadPosRankSum=-4.10000e-01;SOR=1.05100e+00;VQSLOD=-6.57000e+00;VQSR_culprit=MQRankSum;GQ_HIST_ALT=15|19|12|23|42|29|53|93|59|124|167|89|174|172|98|178|139|82|116|8511,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=6|73|421|1185|1785|1676|1300|844|663|514|428|347|263|214|134|120|69|49|35|17,0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|4|9|41|145|387|735|997|1373|1016|1051|570|496|440|365|258|109|31|17|0,0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0;GQ_HIST_ALL=173|42|29|68|96|71|158|218|162|337|412|281|699|353|549|396|680|182|574|9920;DP_HIST_ALL=12|102|557|1623|2757|2714|2369|1963|979|567|442|356|265|216|135|120|69|49|35|18;AB_HIST_ALL=0|4|10|41|145|387|735|995|1373|1018|1052|570|498|440|364|259|109|31|17|0;AC_Male=6665,1;AC_Female=5407,0;AN_Male=16504;AN_Female=13406;AF_Male=4.03841e-01,6.05914e-05;AF_Female=4.03327e-01,0.00000e+00;GC_Male=2756,4321,1172,1,0,0;GC_Female=2189,3619,894,0,0,0;GC_raw=5200,8044,2151,1,0,0;AC_raw=12346,1;AN_raw=30800;GC=4945,7940,2066,1,0,0;AF_raw=4.00844e-01,3.24675e-05;Hom_AFR=631,0;Hom_AMR=106,0;Hom_ASJ=20,0;Hom_EAS=304,0;Hom_FIN=193,0;Hom_NFE=767,0;Hom_OTH=45,0;Hom=2066,0;Hom_raw=2151,0;AC_AFR=3492,0;AC_AMR=442,0;AC_ASJ=91,0;AC_EAS=1019,0;AC_FIN=1460,0;AC_NFE=5207,1;AC_OTH=361,0;AN_AFR=8494;AN_AMR=804;AN_ASJ=260;AN_EAS=1590;AN_FIN=3484;AN_NFE=14330;AN_OTH=948;AF_AFR=4.11114e-01,0.00000e+00;AF_AMR=5.49751e-01,0.00000e+00;AF_ASJ=3.50000e-01,0.00000e+00;AF_EAS=6.40881e-01,0.00000e+00;AF_FIN=4.19059e-01,0.00000e+00;AF_NFE=3.63364e-01,6.97837e-05;AF_OTH=3.80802e-01,0.00000e+00;STAR_AC=4;STAR_AC_raw=5;STAR_Hom=1;POPMAX=EAS,NFE;AC_POPMAX=1019,1;AN_POPMAX=1590,14330;AF_POPMAX=6.40881e-01,6.97837e-05;DP_MEDIAN=30,30;DREF_MEDIAN=5.01187e-55,2.51189e-46;GQ_MEDIAN=99,99;AB_MEDIAN=4.61538e-01,6.00000e-01;AS_RF=7.75613e-01,7.57709e-01;AS_FilterStatus=PASS,PASS;CSQ=A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||1||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001729507|enhancer|||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||,A|intergenic_variant|MODIFIER||||||||||||||||2||||insertion|1||||||||||||||||||||||||||||||||||||||||||||;GC_AFR=1386,2230,631,0,0,0;GC_AMR=66,230,106,0,0,0;GC_ASJ=59,51,20,0,0,0;GC_EAS=80,411,304,0,0,0;GC_FIN=475,1074,193,0,0,0;GC_NFE=2721,3673,767,1,0,0;GC_OTH=158,271,45,0,0,0;Hom_Male=1172,0;Hom_Female=894,0",
            other: b"",
        });
    let mut write_data = Vec::new();
    record.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);
    let complete_data = record.complete_parse().unwrap();
    let mut write_data = Vec::new();
    complete_data.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);

    // 1000genomes
    let data = b"22\t16914962\trs9306245\tC\tA,T\t100\tPASS\tAC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC\tGT\t0|2\t0|0\t0|2\t0|2\t0|2\n";
    let record = PartialVCFRecord::parse_vcf(1, data).unwrap();
    assert_eq!(
        record,
        PartialVCFRecord {
            raw: data,
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 16914962,
            id: Cow::Borrowed(b"rs9306245"),
            reference: Cow::Borrowed(b"C"),
            alternative: vec![Cow::Borrowed(b"A"), Cow::Borrowed(b"T")],
            qual: Cow::Borrowed(b"100"),
            filter: Cow::Borrowed(b"PASS"),
            unparsed_info: Cow::Borrowed(b"AC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC"),
            original_unparsed_info: b"AC=161,1659;AF=0.0321486,0.33127;AN=5008;MULTI_ALLELIC",
            other: b"GT\t0|2\t0|0\t0|2\t0|2\t0|2",
        }
    );
    let mut write_data = Vec::new();
    record.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);
    let complete_data = record.complete_parse().unwrap();
    assert_eq!(
        complete_data,
        CompleteVCFRecord {
            line: 1,
            contig: Cow::Borrowed(b"22"),
            position: 16914962,
            id: Cow::Borrowed(b"rs9306245"),
            reference: Cow::Borrowed(b"C"),
            alternative: vec![Cow::Borrowed(b"A"), Cow::Borrowed(b"T")],
            qual: Cow::Borrowed(b"100"),
            filter: Cow::Borrowed(b"PASS"),
            info: vec![
                (
                    Cow::Borrowed(b"AC"),
                    vec![Cow::Borrowed(b"161"), Cow::Borrowed(b"1659")]
                ),
                (
                    Cow::Borrowed(b"AF"),
                    vec![Cow::Borrowed(b"0.0321486"), Cow::Borrowed(b"0.33127")]
                ),
                (Cow::Borrowed(b"AN"), vec![Cow::Borrowed(b"5008")]),
                (Cow::Borrowed(b"MULTI_ALLELIC"), vec![])
            ],
            format: vec![Cow::Borrowed(b"GT")],
            call: vec![
                vec![vec![Cow::Borrowed(b"0|2")]],
                vec![vec![Cow::Borrowed(b"0|0")]],
                vec![vec![Cow::Borrowed(b"0|2")]],
                vec![vec![Cow::Borrowed(b"0|2")]],
                vec![vec![Cow::Borrowed(b"0|2")]]
            ],
        }
    );
    let mut write_data = Vec::new();
    complete_data.write(&mut write_data).unwrap();
    assert_eq!(&write_data[..], &data[..]);
}

#[test]
fn test_vcf_header_item_parse() {
    let item = VCFHeaderItem::parse(b"##fileformat=VCFv4.1\n", 1).unwrap();
    assert_eq!(
        item,
        VCFHeaderItem {
            key: b"fileformat".to_vec(),
            value: b"VCFv4.1".to_vec(),
            detail: HashMap::new(),
        }
    );

    let mut detail = HashMap::new();
    detail.insert(b"ID".to_vec(), b"1".to_vec());
    detail.insert(b"assembly".to_vec(), b"b37".to_vec());
    detail.insert(b"length".to_vec(), b"249250621".to_vec());
    assert_eq!(
        VCFHeaderItem::parse(b"##contig=<ID=1,assembly=b37,length=249250621>\n", 1).unwrap(),
        VCFHeaderItem {
            key: b"contig".to_vec(),
            value: b"<ID=1,assembly=b37,length=249250621>".to_vec(),
            detail
        }
    );

    let mut detail = HashMap::new();
    detail.insert(b"ID".to_vec(), b"PL".to_vec());
    detail.insert(b"Number".to_vec(), b"G".to_vec());
    detail.insert(b"Type".to_vec(), b"Integer".to_vec());
    detail.insert(
        b"Description".to_vec(),
        b"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"
            .to_vec(),
    );
    assert_eq!(VCFHeaderItem::parse(b"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n", 1).unwrap(), VCFHeaderItem{
        key: b"FORMAT".to_vec(),
        value: b"<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">".to_vec(),
        detail
    });
}

#[test]
fn test_real_gnomad_3_0_merged() -> Result<(), VCFParseError> {
    let mut reader = VCFReader::new(flate2::read::MultiGzDecoder::new(
        &include_bytes!(
            "../../testfiles/gnomad/gnomad.genomes.r3.0.sites.chr22.aconly.subset.merged.vcf.gz"
        )[..],
    ))?;
    let mut counter = 0;
    while let Some(record) = reader.next_record()? {
        record.complete_parse()?;
        counter += 1;
    }
    assert_eq!(counter, 10470);
    Ok(())
}

#[test]
fn test_real_gnomad_3_0_splitted() -> Result<(), VCFParseError> {
    let mut reader = VCFReader::new(flate2::read::MultiGzDecoder::new(
        &include_bytes!(
            "../../testfiles/gnomad/gnomad.genomes.r3.0.sites.chr22.aconly.subset.vcf.gz"
        )[..],
    ))?;
    let mut counter = 0;
    while let Some(record) = reader.next_record()? {
        record.complete_parse()?;
        counter += 1;
    }
    assert_eq!(counter, 12182);
    Ok(())
}

#[test]
fn test_real_1000genomes() -> Result<(), VCFParseError> {
    let mut reader = VCFReader::new(flate2::read::MultiGzDecoder::new(
        &include_bytes!("../../testfiles/1000-genomes/1kGP-subset.vcf.gz")[..],
    ))?;
    let mut counter = 0;
    while let Some(record) = reader.next_record()? {
        record.complete_parse()?;
        counter += 1;
    }
    assert_eq!(counter, 1420);
    Ok(())
}
