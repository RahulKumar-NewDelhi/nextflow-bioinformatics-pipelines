#!/usr/bin/env nextflow

def createSamplesChannel(String filename) {
  def array = []
  def header, line;
  new File(filename).withReader { reader ->
    header = reader.readLine()
    while ((line = reader.readLine()) != null) {
      row = line.split(",")
      def meta = [:]
      meta.id = row[0]
      if (row[1] == "false") {
          meta.pair_end = false
          // optional trailing columns, used only by tumor-normal Mutect2 (SOMATIC_ANALYSIS + params.matchedNormal)
          meta.patient_id = row.size() > 3 ? row[3] : ""
          meta.sample_type = row.size() > 4 ? row[4] : ""
          array << [meta, [file(row[2])]]
      } else if (row[1] == "true"){
          meta.pair_end = true
          meta.patient_id = row.size() > 4 ? row[4] : ""
          meta.sample_type = row.size() > 5 ? row[5] : ""
          array << [meta, [file(row[2]), file(row[3])]]
      }
    }
  }
  return array
}

// Groups sample-sheet rows into tumor/normal pairs for matched-normal Mutect2.
// Requires the sheet's optional trailing patient_id,sample_type columns
// (sample_type values "tumor" / "normal", case-insensitive). Patients missing
// either a tumor or a normal sample are skipped.
def createTumorNormalPairs(String filename) {
  def bySample = [:]
  def header, line;
  new File(filename).withReader { reader ->
    header = reader.readLine()
    while ((line = reader.readLine()) != null) {
      row = line.split(",")
      def sampleId = row[0]
      def paired = row[1]
      def patientIdx = (paired == "true") ? 4 : 3
      def typeIdx = (paired == "true") ? 5 : 4
      def patientId = row.size() > patientIdx ? row[patientIdx].trim() : ""
      def sampleType = row.size() > typeIdx ? row[typeIdx].trim().toLowerCase() : ""
      bySample[sampleId] = [patientId, sampleType]
    }
  }

  def byPatient = [:]
  bySample.each { sampleId, info ->
    def patientId = info[0]
    def sampleType = info[1]
    if (!patientId || !sampleType) return
    byPatient.computeIfAbsent(patientId, { [:] })
    byPatient[patientId][sampleType] = sampleId
  }

  def pairs = []
  byPatient.each { patientId, samples ->
    if (samples.tumor && samples.normal) {
      pairs << [patientId, samples.tumor, samples.normal]
    }
  }
  return pairs
}

def validateSampleSheet(String filename) {
  def array = []
  def header, line;
  new File(filename).withReader { reader ->
    header = reader.readLine()
    while ((line = reader.readLine()) != null) {
      row = line.split(",")
      if (row[1] == "false") {
          fastq = row[2] ==~ /.*.fq|.*.fq.gz|.*.fastq|.*.fastq.gz/
	  array << fastq
      } else if (row[1] == "true"){
          fastq_1 = row[2] ==~ /.*.fq|.*.fq.gz|.*.fastq|.*.fastq.gz/
          fastq_2 = row[3] ==~ /.*.fq|.*.fq.gz|.*.fastq|.*.fastq.gz/
	  array << fastq_1 << fastq_2;
      }
    }
  }
  isValid = (array.isEmpty() == true) ? false: array.every()
  return isValid
}


