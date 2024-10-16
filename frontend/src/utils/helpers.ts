/**
 * The function `readTextFile` reads the text content of a file selected by the user through an HTML
 * input element.
 *
 * Args:
 *   event (Event): The `event` parameter is an event object that is passed to the `readTextFile`
 * function. It is of type `Event` and is used to capture the event that triggered the function, such
 * as a file input change event.
 *
 * Returns:
 *   The `readTextFile` function returns the text content of the file selected by the user in the HTML
 * input element. If a file is selected, the function reads the text content of the file asynchronously
 * and returns it. If no file is selected, an empty string `''` is returned.
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
