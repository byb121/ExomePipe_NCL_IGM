package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.samples.PedigreeValidationType


//import org.broadinstitute.sting.gatk.phonehome._

//import org.broadinstitute.sting.queue.extensions.snpeff._
import org.broadinstitute.sting.queue.util.QScriptUtils

/**
 * HaplotypeCaller
 */
class MyHyplotypeCaller extends QScript {
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="A file contains Bam files to genotype.", shortName="I")
  var bamFile: File = _

  @Input(doc="A file name to for output raw vcf.", shortName="V_out")
  var outVCF: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _


  // This argument is supposed to add padding regions of specified length to the both ends of a splitted region
  // (Queue split the whole genome to the number specified with "sg")
  // but it is not working and GATK developer refused to give support of it, as it is not their implementation.
  @Argument(doc="interval padding ", fullName="interval_padding", shortName="ip",required=false)
  var ip: Int = _

  @Input(doc="An optional file to output active regions.", shortName="activeRegionsOut", required=false)
  var acfile: File = _

  @Argument(fullName="pedigree", shortName="ped", doc="Pedigree files for samples", required=false)
  var ped: List[File] = _
  @Argument(fullName="pedigreeValidationType", shortName="pedValidationType", doc="How strict should we be in validating the pedigree information?",required=false)
  var pedValidationType = PedigreeValidationType.SILENT

  //@Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  //var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.
  //var filterNames: List[String] =  List("LowQualityDepth","MappingQuality","StrandBias","HaplotypeScoreHigh","MQRankSumLow","ReadPosRankSumLow" )

  //@Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  //var filterExpressions: List[String] = Nil
  //var filterExpressions: List[String] = List ("QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0")

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.

  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    // Set the memory limit to 8 gigabytes on each command.
     this.memoryLimit = 8
    //this.javaOpts = ""
    // this.phone_home = GATKRunReport.PhoneHomeOption.NO_ET
    // this.gatk_key = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/TianL_email.chop.edu.key"
    this.interval_padding = ip
    if (qscript.ped == null) {
      this.pedigree =  Nil
    } else {
      this.pedigree =  qscript.ped
      //printf(qscript.ped.toString)
      this.pedValidationType = qscript.pedValidationType
    }
  }

  val queueLogDir: String = ".qlog/" // Gracefully hide Queue's output

  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  def script() {
    // Create the four functions that we may run depending on options.
    val genotyper = new HaplotypeCaller with UnifiedGenotyperArguments
    //val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    //val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    val bams = QScriptUtils.createSeqFromFile(bamFile)

    genotyper.scatterCount = nContigs
    // genotyper.input_file = qscript.bamFile
    genotyper.input_file = bams
    genotyper.out = outVCF
    genotyper.activeRegionOut = acfile

    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")
    evalUnfiltered.num_threads = 1

    //variantFilter.variant = genotyper.out
    //variantFilter.out = swapExt(genotyper.out, "vcf", "filtered.vcf")
    //variantFilter.filterName = filterNames
    //variantFilter.filterExpression = filterExpressions

    //evalFiltered.eval :+= variantFilter.out
    //evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")
    //evalFiltered.num_threads = 1

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    //if (filterNames.size > 0)
    // add(variantFilter, evalFiltered)

    // annotate_snpEff ( variantFilter.out )

  }



  /*case class varannotator (inVcf: File, inSnpEffFile: File, outVcf: File) extends VariantAnnotator  {
      this.variant = inVcf
      this.snpEffFile = inSnpEffFile
      this.out = outVcf
      this.alwaysAppendDbsnpId = true
      this.D = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"
      this.R = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/hg19/hg19.fa"
      this.A = Seq("SnpEff")
      this.isIntermediate = false
      this.analysisName = queueLogDir + outVcf + ".varannotator"
      this.jobName = queueLogDir + outVcf + ".varannotator"
      this.scatterCount = nContigs
    }


    def annotate_snpEff(inVcf: File) {
          val eff = new SnpEff
          eff.config = new File("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.config")
          eff.genomeVersion = "GRCh37.64"

          eff.inVcf = inVcf
          var snpEffout: File  = swapExt(eff.inVcf, "vcf", "snpEff.vcf")
          eff.outVcf = swapExt(eff.inVcf, "vcf", "snpEff.out")

          eff.javaClasspath = List("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/")
          eff.jarFile = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.jar"

          add(eff)
          add(varannotator(eff.inVcf,eff.outVcf,snpEffout))
  }*/


}



