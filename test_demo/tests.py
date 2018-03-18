from django.test import TestCase
from .models import Animal

# Create your tests here.

class AnimalTestCase(TestCase):
    def setUp(self):
        Animal.objects.create(name="lion",sound="roar")
        Animal.objects.create(name="cat", sound="meow")
    def test_animals_can_speak(self):
        """Animals that can speak are correctly identified"""
        lion = Animal.objects.get(name="lion")
        cat = Animal.objects.get(name="cat")
        self.assertEqual(lion.speak(), 'The lion says "roar"')
        self.assertEqual(cat.speak(), 'The cat says "meow"')
    #./manage.py test test_demo.tests.AnimalTestCase
