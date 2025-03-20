export function checkFasta(content: string) {
  const errors = []

  if (content == '') {
    errors.push('File is empty.')
  } else {
    const sequenceBlocks = content.trim().split(/\n{2,}/)

    for (const block of sequenceBlocks) {
      const parsedBlock = block.match(/^\>(\S+) ?.*?\r?\n(.+)/s)

      const header = parsedBlock?.[1]
      const sequence = parsedBlock?.[2]

      if (header && sequence) {
        const invalidCharacters = sequence?.match(/[^ACTGU\r\n]/g)

        if (invalidCharacters) {
          for (const match of invalidCharacters) {
            errors.push(`Invalid character ${match[0]} for ${header}.`)
          }
        }

        const flatSequence = sequence.replace(/\s/g, '')

        if (flatSequence.length % 3 !== 0) {
          errors.push(
            `Invalid number of nucleotides (${flatSequence.length}) for ${header}. Should be a multiple of three.`,
          )
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

export function checkClustal(content: string) {
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
) {
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

export function checkUtrSequence(sequence: string) {
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
          errors.push(`Invalid character ${matchedChar}.`)
          break
      }
    }
  }

  return errors
}

export function checkPasswordConstraints(password: string) {
  if (password.length < 8 || password.length > 20) {
    return 'The password must be between 8 and 20 characters long.'
  } else if (!password.match(/([a-z].*\d)|(\d.*[a-z])/i)) {
    return 'The password must contain at least one letter and one number.'
  } else {
    return ''
  }
}
