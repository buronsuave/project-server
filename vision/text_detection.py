import os
import base64

from google.cloud import vision

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = r'ServiceAccount.json'
client = vision.ImageAnnotatorClient()

def detectText(imgText):
    imgTextEnc = imgText.encode('utf-8')
    decoded_image_data = base64.decodebytes(imgTextEnc)

    image = vision.Image(content = decoded_image_data)
    textDetected = client.document_text_detection(image = image)

    finalText = textDetected.full_text_annotation.text

    try:
        endRow = finalText.index('\n')
    except:
        endRow = -1
    
    if endRow != -1:
        finalText = finalText[:endRow]
    
    finalText = finalText.replace(" ", "")
    finalText = finalText.replace("Dl", "D(")
    return finalText

