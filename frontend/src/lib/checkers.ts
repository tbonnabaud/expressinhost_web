export function checkFasta(content: string) {
  const errors = []

  if (content == '') {
    errors.push('File is empty')
  }

  return errors
}

export function checkClustal(content: string) {
  const errors = []

  if (content == '') {
    errors.push('File is empty')
  }

  return errors
}
