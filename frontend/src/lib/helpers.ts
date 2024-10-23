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
