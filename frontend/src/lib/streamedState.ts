import { ref, onUnmounted } from 'vue'

export enum Status {
  IDLE = 'Idle',
  RUNNING = 'Running',
  ERROR = 'Error',
  SUCCESS = 'Success',
}

export interface StreamState {
  status: Status
  message: string
  done: number | null
  total: number | null
  result?: object | null
}

export function useStreamState(url: string, method: string, token?: string) {
  const state = ref<StreamState | null>(null)
  const isStreaming = ref(false)
  let controller: AbortController | null = null

  const startStream = async (data?: object) => {
    isStreaming.value = true
    controller = new AbortController()
    const signal = controller.signal
    const headers = {
      Accept: 'application/json',
      'Content-Type': 'application/json',
    }

    try {
      const response = await fetch(url, {
        method: method,
        headers: token
          ? { ...headers, Authorization: `Bearer ${token}` }
          : headers,
        body: data ? JSON.stringify(data) : undefined,
        signal,
      })

      if (!response.body) throw new Error('Stream not available')

      const reader = response.body.getReader()
      const decoder = new TextDecoder()
      let partialData = ''

      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        partialData += decoder.decode(value, { stream: true })

        try {
          // Try to parse data
          const parsed: StreamState = JSON.parse(partialData)
          state.value = parsed
          // Reset partial data after successful parsing
          partialData = ''
        } catch (parseError) {
          console.error(parseError)
        }
      }
    } catch (error) {
      if (error instanceof Error) {
        if (error.name === 'AbortError') {
          console.log('Stream aborted')
        } else {
          console.error('Error fetching data:', error.message)
        }
      }
    } finally {
      isStreaming.value = false
    }
  }

  const stopStream = () => controller?.abort()

  // onMounted(() => {
  //   startStream()
  // })

  onUnmounted(() => {
    stopStream()
  })

  return { state, isStreaming, startStream, stopStream }
}
