// Stop codons for DNA and mRNA
const STOP_CODONS = ['UAG', 'UAA', 'UGA', 'TAG', 'TAA', 'TGA']

export function checkSequence(sequence: string): string[] {
  const errors: string[] = []

  const invalidCharacters = sequence?.match(/[^ACTGU\r\n]/g)

  if (invalidCharacters) {
    for (const match of invalidCharacters) {
      errors.push(`Invalid character "${match[0]}".`)
    }
  }

  // Flatten the sequence
  const flatSequence = sequence.replace(/\s/g, '')

  if (flatSequence.length % 3 !== 0) {
    errors.push(
      `Invalid number of nucleotides (${flatSequence.length}). Should be a multiple of three.`,
    )
  } else {
    // Check if there are stop codons before the end of the sequence
    for (let i = 0; i < flatSequence.length - 3; i += 3) {
      const codon = flatSequence.substring(i, i + 3)

      if (STOP_CODONS.includes(codon)) {
        errors.push(`Stop codon "${codon}" before the end of the sequence.`)
      }
    }
  }

  return errors
}

export function checkFasta(content: string): string[] {
  const errors: string[] = []

  if (content == '') {
    errors.push('File is empty.')
  } else {
    const sequenceBlocks = content.trim().split(/\n{2,}/)

    for (const block of sequenceBlocks) {
      const parsedBlock = block.match(/^\>(\S+) ?.*?\r?\n(.+)/s)

      const header = parsedBlock?.[1]
      const sequence = parsedBlock?.[2]

      if (header && sequence) {
        const sequenceErrors = checkSequence(sequence)

        for (const seqError of sequenceErrors) {
          errors.push(`${header}: ${seqError}`)
        }
      } else if (header && !sequence) {
        errors.push(`Missing sequence for ${header}.`)
      } else {
        errors.push('Invalid format.')
        // If one block is not correctly formatted return directly the result
        return errors
      }
    }
  }

  return errors
}

/**
 * Function to validate if a given content is a correct PDB format.
 * @param content - The content to validate as a PDB.
 * @returns A list of string errors indicating issues with the PDB format.
 */
export function checkPdb(content: string): string[] {
  const errors: string[] = []
  const lines = content.split('\n')

  // Check if the file starts with the correct header
  if (!lines[0]?.startsWith('HEADER')) {
    errors.push('Missing or incorrect HEADER record.')
  }

  // Check for required records in the PDB file
  const requiredRecords = ['TITLE', 'COMPND', 'SOURCE', 'SEQRES', 'ATOM', 'END']

  for (const record of requiredRecords) {
    if (!lines.some(line => line.startsWith(record))) {
      errors.push(`Missing required record: ${record}`)
    }
  }

  // Additional validation can be added here, such as checking the format of each line
  // For simplicity, we'll assume the presence of required records is sufficient

  return errors
}

export function checkClustal(content: string): string[] {
  const errors = []

  if (content == '') {
    errors.push('File is empty.')
  } else {
    if (!content.match(/^CLUSTAL/)) {
      errors.push('Missing "CLUSTAL" keyword at the beginning of the file.')
    }
  }

  return errors
}

export function checkClustalMatchingFasta(
  clustalContent: string,
  fastaContent: string,
): string[] {
  const errors = []
  const fastaIDs = Array.from(fastaContent.matchAll(/^>(\S+)/gm), m => m[1])
  const clustalBlocks = clustalContent
    .trim()
    .split(/\n{2,}/)
    .slice(1) // Slice to ignore the block with "CLUSTAL" header

  for (const [blockIndex, block] of clustalBlocks.entries()) {
    const clustalLines = Array.from(
      block.matchAll(/^(\S+)\s+([A-Z\-]+)\s+(\d+)/gm),
    )
    const blockNumber = blockIndex + 1

    if (clustalLines.length !== fastaIDs.length) {
      errors.push(
        `Wrong number of lines in block ${blockNumber}, should be ${fastaIDs.length}.`,
      )
    } else {
      for (const [lineIndex, line] of clustalLines.entries()) {
        const clustalID = line[1]

        if (fastaIDs[lineIndex] !== clustalID) {
          errors.push(
            `IDs matching error in block ${blockNumber}. FASTA: ${fastaIDs[lineIndex]}, CLUSTAL: ${clustalID}.`,
          )
        }
      }
    }
  }

  return errors
}

export function checkUtrSequence(sequence: string): string[] {
  const errors = []
  const trimmedSequence = sequence.trim()

  const invalidCharacters = trimmedSequence.match(/[^ACTGU]/g)

  if (invalidCharacters) {
    for (const match of invalidCharacters) {
      const matchedChar = match[0]

      switch (matchedChar) {
        case '\r\n':
        case '\n':
          errors.push(`Invalid newline character.`)
          break
        case ' ':
          errors.push(`Invalid space character.`)
          break
        default:
          errors.push(`Invalid character "${matchedChar}".`)
          break
      }
    }
  }

  return errors
}

export function checkPasswordConstraints(password: string): string {
  if (password.length < 8 || password.length > 20) {
    return 'The password must be between 8 and 20 characters long.'
  } else if (!password.match(/([a-z].*\d)|(\d.*[a-z])/i)) {
    return 'The password must contain at least one letter and one number.'
  } else {
    return ''
  }
}
