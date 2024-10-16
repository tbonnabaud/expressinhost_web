/**
 * This TypeScript function reads the text content of a file selected by the user through an HTML input
 * element.
 * @param {Event} event - The `event` parameter in the `readTextFile` function is an event object that
 * is passed to the function when it is triggered. In this case, it is expected to be of type `Event`,
 * which is a standard DOM event object. The event object contains information about the event that
 * occurred.
 * @returns The `readTextFile` function returns the text content of the file selected by the user if a
 * file is selected. If no file is selected, it returns an empty string.
 */
export async function readTextFile(event: Event) {
  const target = event.target as HTMLInputElement

  if (target.files && target.files.length > 0) {
    const file = target.files[0]
    const text = await file.text()
    return text
  }

  return ''
}
