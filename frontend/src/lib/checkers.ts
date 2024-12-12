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
