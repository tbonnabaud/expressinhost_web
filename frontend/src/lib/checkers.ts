export function checkFasta(content: string) {
  const errors = []

  if (content == '') {
    errors.push('File is empty.')
  } else {
    const sequenceBlocks = content.trim().split(/\n{2,}/)

    for (const block of sequenceBlocks) {
      const parsedBlock = block.match(/^\>(\S+) ?.*?\n(.+)/s)

      const header = parsedBlock?.[1]
      const sequence = parsedBlock?.[2]

      if (header && sequence) {
        const invalidCharacters = sequence?.match(/[^ACTGU\n]/g)

        if (invalidCharacters) {
          for (const match of invalidCharacters) {
            errors.push(`Invalid character ${match[0]} for ${header}.`)
          }
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
