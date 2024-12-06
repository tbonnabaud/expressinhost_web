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
