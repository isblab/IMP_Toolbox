```python
def make_protein(protein_name, seq_start, seq_end, *kwargs):
    template_dict = {
            "mode": "protein",
            "config": {
            "canvasWidth": seq_end - seq_start + 1,
            "canvasHeight": 600,
            "viewBox": "10 -230 1000 600",
            "scaleLength": 3000,
            "zoomRatio": 1
            },
            "data":[
                {
                "type": "protein",
                "length": seq_end - seq_start + 1,
                "option": {
                    "displayName": protein_name,
                    "position": {
                        "start": {
                            "site": seq_start,
                            "display": "top"
                        },
                        "end": {
                            "site": seq_end,
                            "display": "top"
                        }
                    },
                    "coordinate": {
                        "vertical": {
                            "start": 1,
                            "isLocked": False
                        },
                        "horizontal": {
                            "start": 0,
                            "end": seq_end - seq_start,
                            "isLocked": False
                        }
                    },
                    "id": "Protein_" + protein_name,
                    "style":{
                        "align": "custom",
                        "height": 25,
                        "fontSize": 12,
                        "color": "#babdb6",
                        "gradient": "none",
                        "texture": {
                        "type": "none",
                        "color": "#333333"
                        }
                    },
                    "borderStyle": {
                        "color": "#000000",
                        "size": 1,
                        "isDash": False
                    }
                },
                "children": []
            }
            ]
     }

    return template_dict
```