import type { CodonTranslation } from './interfaces'

/**
 * The function `readTextFile` reads the text content of a file selected by the user through an HTML
 * input element.
 * @param {Event} event - Event object that represents the event that triggered the function.
 * @returns The function `readTextFile` returns the text content of the file selected by the user if a
 * file is selected. If no file is selected, it logs an error message and returns an empty string.
 */
export async function readTextFile(event: Event) {
  const target = event.target as HTMLInputElement

  if (target.files && target.files.length > 0) {
    const file = target.files[0]
    const text = await file.text()

    return text
  }

  console.error('No file.')

  return ''
}

/**
 * The function `sleep` returns a promise that resolves after a specified number of milliseconds.
 * @param {number} ms - The `ms` parameter in the `sleep` function represents the number of
 * milliseconds for which the function will pause execution before resolving the promise.
 * @returns The `sleep` function is returning a Promise.
 */
export function sleep(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms))
}

/**
 * The function `formatToLocaleDateString` converts a UTC date string to a localized date string.
 * @param {string} utcDateString - A string representing a date and time in Coordinated Universal Time
 * (UTC) format.
 * @returns The function `formatToLocaleDateString` takes a UTC date string as input, converts it to a
 * Date object, and then returns the date in a localized string format using the `toLocaleString`
 * method.
 */
export function formatToLocaleDateString(utcDateString: string) {
  return new Date(utcDateString).toLocaleString()
}

/**
 * The `groupByAminoAcid` function takes an array of `CodonTranslation` objects and groups them by
 * their `amino_acid` property.
 * @param {CodonTranslation[]} array - The `array` parameter in the `groupByAminoAcid` function is an
 * array of `CodonTranslation` objects. Each `CodonTranslation` object likely represents a codon and
 * its corresponding amino acid.
 * @returns The `groupByAminoAcid` function is returning an object where the keys are amino acids and
 * the values are arrays of `CodonTranslation` objects that correspond to that amino acid.
 */
export function groupByAminoAcid(array: CodonTranslation[]) {
  return array.reduce(
    (acc, e) => {
      const aminoAcid = e.amino_acid
      acc[aminoAcid] = acc[aminoAcid] ? [...acc[aminoAcid], e] : [e]
      return acc
    },
    {} as Record<string, CodonTranslation[]>,
  )
}

/**
 * The `toFixedFloat` function rounds a number to a specified number of decimal places
 * and returns it as a floating-point number.
 */
export function toFixedFloat(value: number, fractionDigits: number = 0) {
  return parseFloat(value.toFixed(fractionDigits))
}

/**
 * The `downloadFile` function in TypeScript creates a downloadable file from the provided content and
 * filename.
 * @param {string} content - The `content` parameter in the `downloadFile` function is the actual data
 * that you want to download as a file. It should be a string containing the content that you want to
 * save in the file.
 * @param {string} filename - The `filename` parameter is a string that represents the name under which
 * the file will be saved when downloaded. For example, if you are downloading a text file, the
 * `filename` could be something like "example.txt".
 */
export function downloadFile(content: string, filename: string) {
  // Creating a Blob from the data
  const blob = new Blob([content], {
    type: 'text/plain',
  })
  const url = URL.createObjectURL(blob)

  // Creating a temporary link element
  const link = document.createElement('a')
  link.href = url
  link.download = filename

  // Append to the body, click and remove it
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)

  // Revoke the object URL after the download
  URL.revokeObjectURL(url)
}

/**
 * Parses an ISO 8601 duration string and returns the total duration in seconds.
 *
 * @param duration - The ISO 8601 duration string to parse.
 *   The format should be "PT[nH][nM][nS]", where:
 *   - "PT" is the time designator for time durations.
 *   - "nH" represents hours (optional).
 *   - "nM" represents minutes (optional).
 *   - "nS" represents seconds (optional).
 *   Examples: "PT1.5S", "PT1H30M15.5S".
 *
 * @returns The total duration in seconds.
 *
 * @example
 * parseISODuration("PT1.565351S"); // returns 1.565351
 * parseISODuration("PT1H30M15.5S"); // returns 5415.5
 */
export function parseISODuration(duration: string): number {
  const regex = /(?:(\d+)H)?(?:(\d+)M)?(?:(\d+(\.\d+)?)S)?/
  const matches = duration.match(regex)

  if (!matches) {
    return 0
  }

  const hours = parseInt(matches[1] || '0')
  const minutes = parseInt(matches[2] || '0')
  const seconds = parseFloat(matches[3] || '0')

  const totalSeconds = seconds + minutes * 60 + hours * 3600

  return totalSeconds
}
